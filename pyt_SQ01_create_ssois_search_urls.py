import re
import sys
import time
import datetime
import math
import numpy as np
import urllib
import glob, os, bz2, subprocess
import os.path
from astropy.time import Time, TimeDelta
from astropy.io.votable import parse_single_table, VOTableSpecWarning, VOWarning
from decimal import *
import sqlite3
from sqlite3 import Error
import jdcal
from jdcal import gcal2jd,jd2gcal
import requests
import warnings
import statistics
from uncertainties import unumpy,ufloat
from uncertainties.umath import *
from uncertainties.unumpy import *
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- URL GENERATION FUNCTIONS #####

def retrieve_instrument_id(instrument_name,sqlite_file,keep_going,path_logfile,path_errorfile):
    # Retrieve instrument_id corresponding to a given instrument name
    instrument_id = int(0)
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving instrument ID for {:s} from database...'.format(instrument_name))
        try:
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            query = "SELECT instrument_id FROM instruments WHERE instrument_name='{:s}'".format(instrument_name)
            mcd.output_log_entry(path_logfile,query)
            cursor.execute(query)
            row = cursor.fetchone()
            if row != None:
                instrument_id = int(row[0])
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Instrument {:s} not found'.format(instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for instrument {:s}: retrieve_instrument_id()'.format(instrument_name))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_instrument_id()')
    return instrument_id,keep_going


def generate_search_urls_numbrasts(sqlite_file,instrument_id,current_datetime,search_urls,search_history_file,telinst,keep_going,path_logfile,path_errorfile):
    # Generate search URLS for numbered asteroids
    if keep_going:
        try:
            if telinst == 'SDSS':    telinst_search_term = 'SDSS'
            if telinst == 'MegaCam': telinst_search_term = 'CFHT%2FMegaCam'
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            mcd.output_log_entry(path_logfile,'Generating search URLs for numbered asteroids...')
            mcd.output_log_entry(path_logfile,'Retrieving list of numbered asteroids for searching...')
            query = "SELECT sso.ssobject_id,sso.desig_number,sh.last_searched FROM ssobjects AS sso LEFT OUTER JOIN search_history AS sh ON sso.ssobject_id=sh.ssobject_id WHERE sso.desig_number IS NOT NULL AND sso.object_type='asteroid' AND (sh.instrument_id={:d} OR sh.instrument_id IS NULL)".format(instrument_id)
            mcd.output_log_entry(path_logfile,query)
            cursor.execute(query)
            rows = cursor.fetchall()
            for row in rows:
                ssobject_id   = row[0]
                desig_number  = row[1].lstrip()
                last_searched = row[2]
                epoch1 = '1990+01+01'
                epoch2 = '{:4d}+{:02d}+{:02d}'.format(int(current_datetime[0:4]),int(current_datetime[4:6]),int(current_datetime[6:8]))
                current_datetime1 = datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
                mjd_current_date  = gcal2jd(int(current_datetime1[0:4]),int(current_datetime1[4:6]),int(current_datetime1[6:8]))[1]
                if last_searched != None:
                    mjd_searched = gcal2jd(int(last_searched[0:4]),int(last_searched[5:7]),int(last_searched[8:10]))[1]
                    days_since_searched = int(mjd_current_date - mjd_searched)
                    if days_since_searched < 365:
                        mcd.output_log_entry(path_logfile,'Search conducted for asteroid {:s} in {:s} data within last year ({:d} days ago)'.format(desig_number,telinst,days_since_searched))
                        search_desired = False
                    else:
                        # if prior entry for object/instrument combination is >1 year ago, conduct new search
                        mcd.output_log_entry(path_logfile,'Search not conducted for asteroid {:s} in {:s} data within last year (last search {:d} days ago)'.format(desig_number,telinst,days_since_searched))
                        mjd_epoch1  = mjd_searched + 1
                        gcal_epoch1 = jd2gcal(2400000.5,mjd_epoch1)
                        epoch1 = '{:4d}+{:02d}+{:02d}'.format(int(gcal_epoch1[0]),int(gcal_epoch1[1]),int(gcal_epoch1[2]))
                        search_desired = True
                else:
                    # if no prior entry for object/instrument combination is found, conduct new search
                    search_desired = True
                if search_desired:
                    search_urls.write('http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en&object={:d}&search=bynameCADC&epoch1={:s}&epoch2={:s}&eellipse=&eunits=arcseconds&extres=yes&xyres=yes&telinst={:s}\n'.format(int(desig_number),epoch1,epoch2,telinst))
                    search_history_file.write('{:>12d}   {:>13d}   {:s}\n'.format(ssobject_id,instrument_id,datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            conn.close()  # Close connection to database file
            mcd.output_log_entry(path_logfile,'SSOIS search URLs for numbered asteroids written to file')
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: generate_search_urls_numbrasts()')
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: generate_search_urls_numbrasts()')
    return search_urls,keep_going


def split_search_urls(search_urls_master_filename,num_threads,keep_going,path_logfile,path_errorfile):
    search_urls_files = ['' for idx in range(num_threads)]
    for idx in range(0,num_threads):
        search_urls_files[idx] = search_urls_master_filename[:-11] + '_{:02d}.txt'.format(idx+1)
    num_urls = sum(1 for line in open(search_urls_master_filename)) - 1
    num_urls_thread = int(num_urls/num_threads) + 1
    
    idx_output_file,idx_url_in_current_file = 0,0
    with open(search_urls_master_filename,'r') as search_urls_master_file:
        #for _ in range(1): #skip first 2 header lines
        #    next(search_urls_file)
        line = search_urls_master_file.readline()  # copy first line of master URL file to all individual URL files
        for idx in range(0,num_threads):
            with open(search_urls_files[idx],'w') as of:
                of.write(line[:-2] + ' - Part {:d} of {:d}\n'.format(idx+1,num_threads))
        for line in search_urls_master_file:
            with open(search_urls_files[idx_output_file],'a') as of:
                of.write(line)
            idx_url_in_current_file += 1
            if idx_url_in_current_file == num_urls_thread:
                idx_output_file += 1
                idx_url_in_current_file = 0
    
    return search_urls_files,keep_going


def main():
    # Define filenames and paths
    if len(sys.argv)!=5:
        print('Usage:\n python3 pyt_SQ01_create_ssois_search_urls.py [base_path] [sqlite_file] [telinst] [num_threads]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    telinst     = sys.argv[3]
    num_threads = int(sys.argv[4])

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('SQ01_create_ssois_search_urls execution started','SQ01_create_ssois_search_urls execution started.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'SQ01_create_ssois_search_urls')

        # Connect to database file
        ssois_directory = base_path + 'ssois_queries/'
        search_urls_directory = ssois_directory + 'search_urls/'
        mcd.create_directory(ssois_directory,path_logfile,path_errorfile)
        mcd.create_directory(search_urls_directory,path_logfile,path_errorfile)

        #num_threads = 10
        instrument_name = ''
        if telinst == 'SDSS': instrument_name = 'SDSS Imaging Camera'
        if telinst == 'MegaCam': instrument_name = 'MegaCam'
        current_datetime1 = datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        search_urls_master_filename      = search_urls_directory + 'search_urls_{:s}_{:s}_master.txt'.format(telinst,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        search_history_toingest_filename = search_urls_directory + 'search_history_{:s}_{:s}_toingest.txt'.format(telinst,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    
        os.chdir(search_urls_directory)
        with open(search_urls_master_filename,'w') as search_urls_master_file, open(search_history_toingest_filename,'w') as search_history_file:
            search_urlfiles_datetime = datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')
            search_urls_master_file.write('SSOIS search URLs for {:s}\n'.format(search_urlfiles_datetime))
            search_history_file.write(' ssobject_id   instrument_id   last_searched\n'.format(search_urlfiles_datetime))
            instrument_id,keep_going = retrieve_instrument_id(instrument_name,sqlite_file,keep_going,path_logfile,path_errorfile)
            search_urls,keep_going   = generate_search_urls_numbrasts(sqlite_file,instrument_id,current_datetime1,search_urls_master_file,search_history_file,telinst,keep_going,path_logfile,path_errorfile)

        search_urls_files,keep_going = split_search_urls(search_urls_master_filename,num_threads,keep_going,path_logfile,path_errorfile)
        print('Search URLs written to:')
        for idx in range(0,num_threads):
            print(' {:s}'.format(search_urls_files[idx]))
            
        mcd.compress_file_gzip(search_urls_master_filename)
        for idx in range(0,num_threads):
            mcd.compress_file_gzip(search_urls_files[idx])
        mcd.compress_file_gzip(search_history_toingest_filename)
            
    mcd.send_status_email('SQ01_create_ssois_search_urls execution complete','SQ01_create_ssois_search_urls execution complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None


if __name__ == '__main__':
    main()

