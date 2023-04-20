import re
import sys
import time
import datetime
import math
import numpy as np
import urllib.request
import glob, os, bz2, subprocess
import os.path
from numpy import linspace
from decimal import *
import sqlite3
from sqlite3 import Error
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- ORBITAL ELEMENT FILE PROCESSING #####

def download_jpl_orbelem_files(numbr_url,unnum_url,comet_url,numbr_filepath,unnum_filepath,comet_filepath,path_logfile,path_errorfile):
    # Download orbital element files
    mcd.output_log_entry(path_logfile,'Downloading orbital element files...')
    urllib.request.urlretrieve(numbr_url,numbr_filepath)
    urllib.request.urlretrieve(unnum_url,unnum_filepath)
    urllib.request.urlretrieve(comet_url,comet_filepath)
    mcd.output_log_entry(path_logfile,'Downloads complete')
    return None
    
def ingest_numbr_data(numbr_filepath,sqlite_file,path_logfile,path_errorfile):
    # Ingesting orbital element data for numbered asteroids
    mcd.output_log_entry(path_logfile,'Ingesting orbital element data for numbered asteroids from {:s}...'.format(numbr_filepath))
    # Connect to database file
    conn = mcd.create_connection(sqlite_file)
    cursor = conn.cursor()
    with open(numbr_filepath) as input_file:
        for _ in range(2): #skip first two header lines
            next(input_file)
        for line in input_file:
            asteroid_number = int(line[0:6].lstrip())
            asteroid_name = mcd.add_apostrophe_escape(line[7:24].rstrip())
            epoch       = int(line[25:30])
            orbelems    = line[31:].split()
            semimaj     = float(orbelems[0])
            eccen       = float(orbelems[1])
            inclin      = float(orbelems[2])
            argperi     = float(orbelems[3])
            longascnode = float(orbelems[4])
            meananom    = float(orbelems[5])
            hmag        = float(orbelems[6])
            gparam      = float(orbelems[7])
            
            mcd.output_log_entry(path_logfile,'Adding/updating data for asteroid ({:d}) {:s}...'.format(asteroid_number,asteroid_name))
            
            # Check if unnumbered asteroid entry exists with same provisional designation,
            #  and if so, add numerical designation
            query = "SELECT ssobject_id,desig_number,desig_provisional FROM ssobjects WHERE desig_provisional='{:>10s}' and object_type='asteroid'".format(asteroid_name)
            cursor.execute(query)
            mcd.output_log_entry(path_logfile,query)
            row = cursor.fetchone()
            if row != None and row[1] == None:
                query = "UPDATE ssobjects SET desig_number='{:>7d}' WHERE ssobject_id={:d}".format(asteroid_number,row[0])
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            
            # Insert new asteroid entries from newly downloaded JPL catalog
            if asteroid_name[0].isalpha():
                query = "INSERT OR IGNORE INTO ssobjects(object_type,desig_number,desig_name,h_mag,g_param,epoch_mjd,semimaj_axis,eccentricity,inclination,arg_perihelion,long_ascnode,mean_anomaly,date_added) VALUES ('asteroid','{:>7d}','{:>17s}',{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},'{:s}')".format(asteroid_number,
                    asteroid_name,hmag,gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            else:
                query = "INSERT OR IGNORE INTO ssobjects(object_type,desig_number,desig_provisional,h_mag,g_param,epoch_mjd,semimaj_axis,eccentricity,inclination,arg_perihelion,long_ascnode,mean_anomaly,date_added) VALUES ('asteroid','{:>7d}','{:>10s}',{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},'{:s}')".format(asteroid_number,
                    asteroid_name,hmag,gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            
            # Update existing asteroid entries from newly downloaded JPL catalog
            if asteroid_name[0].isalpha():
                query = "UPDATE ssobjects SET object_type='asteroid',desig_name='{:>17s}',h_mag={:f},g_param={:f},epoch_mjd={:f},semimaj_axis={:f},eccentricity={:f},inclination={:f},arg_perihelion={:f},long_ascnode={:f},mean_anomaly={:f},date_updated='{:s}' WHERE desig_number='{:>7d}'".format(asteroid_name,
                    hmag,gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),asteroid_number)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            else:
                query = "UPDATE ssobjects SET object_type='asteroid',desig_provisional='{:>10s}',h_mag={:f},g_param={:f},epoch_mjd={:f},semimaj_axis={:f},eccentricity={:f},inclination={:f},arg_perihelion={:f},long_ascnode={:f},mean_anomaly={:f},date_updated='{:s}' WHERE desig_number='{:>7d}'".format(asteroid_name,
                    hmag,gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),asteroid_number)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)

    # Committing changes and closing the connection to the database file
    conn.commit()
    conn.close()
    mcd.output_log_entry(path_logfile,'Insertion of osculating orbital elements for numbered asteroids from {:s} complete'.format(numbr_filepath))
    return None


def ingest_unnum_data(unnum_filepath,sqlite_file,path_logfile,path_errorfile):
    # Ingesting orbital element data for unnumbered asteroids
    mcd.output_log_entry(path_logfile,'Ingesting orbital element data for unnumbered asteroids from {:s}...'.format(unnum_filepath))
    # Connect to database file
    conn = mcd.create_connection(sqlite_file)
    cursor = conn.cursor()
    with open(unnum_filepath) as input_file:
        for _ in range(2): #skip first two header lines
            next(input_file)
        for line in input_file:
            asteroid_desig = mcd.add_apostrophe_escape(line[0:10].rstrip())
            epoch       = int(line[12:17])
            orbelems    = line[18:].split()
            semimaj     = float(orbelems[0])
            eccen       = float(orbelems[1])
            inclin      = float(orbelems[2])
            argperi     = float(orbelems[3])
            longascnode = float(orbelems[4])
            meananom    = float(orbelems[5])
            hmag        = float(orbelems[6])
            gparam      = float(orbelems[7])

            mcd.output_log_entry(path_logfile,'Adding/updating data for asteroid {:s}...'.format(asteroid_desig))
            # Check if asteroid entry exists with same provisional designation
            query = "SELECT ssobject_id FROM ssobjects WHERE desig_provisional='{:>10s}' and object_type='asteroid'".format(asteroid_desig)
            mcd.output_log_entry(path_logfile,query)
            cursor.execute(query)
            row = cursor.fetchone()
            if row == None:  # If no asteroid entry with same provisional designation exists, add new entry
                query = "INSERT OR IGNORE INTO ssobjects(object_type,desig_provisional,h_mag,g_param,epoch_mjd,semimaj_axis,eccentricity,inclination,arg_perihelion,long_ascnode,mean_anomaly,date_added) VALUES ('asteroid','{:>10s}',{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},'{:s}')".format(asteroid_desig,
                    hmag,gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            else:  # If asteroid entry with same provisional designation exists, update entry
                query = "UPDATE ssobjects SET object_type='asteroid',h_mag={:f},g_param={:f},epoch_mjd={:f},semimaj_axis={:f},eccentricity={:f},inclination={:f},arg_perihelion={:f},long_ascnode={:f},mean_anomaly={:f},date_updated='{:s}' WHERE desig_provisional='{:>10s}'".format(hmag,
                    gparam,epoch,semimaj,eccen,inclin,argperi,longascnode,meananom,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),asteroid_desig)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
    # Committing changes and closing the connection to the database file
    conn.commit()
    conn.close()
    mcd.output_log_entry(path_logfile,'Insertion of osculating orbital elements for unnumbered asteroids from {:s} complete.'.format(unnum_filepath))
    return None


def ingest_comet_data(comet_filepath,sqlite_file,path_logfile,path_errorfile):
    # Ingesting orbital element data for comets
    mcd.output_log_entry(path_logfile,'Ingesting orbital element data for comets from {:s}...'.format(comet_filepath))
    # Connect to database file
    conn = mcd.create_connection(sqlite_file)
    cursor = conn.cursor()
    with open(comet_filepath) as input_file:
        for _ in range(2): #skip first two header lines
            next(input_file)
        for line in input_file:
            if line[0:3].lstrip() != '':
                comet_number = line[0:4]
                comet_desig  = mcd.add_apostrophe_escape(line[5:43].rstrip())
                # Add letter designations for numbered comets with multiple components
                if comet_desig[-3] == '-':
                    comet_number = comet_number + comet_desig[-3:]
                elif comet_desig[-2] == '-':
                    comet_number = comet_number + comet_desig[-2:]
            else:
                comet_number = None
                comet_desig  = mcd.add_apostrophe_escape(line[3:43].rstrip())
            epoch        = int(line[44:51])
            orbelems     = line[52:].split()
            peridist     = float(orbelems[0])
            eccen        = float(orbelems[1])
            inclin       = float(orbelems[2])
            argperi      = float(orbelems[3])
            longascnode  = float(orbelems[4])
            time_peri    = float(orbelems[5])

            if comet_number != None:  # If input comet is numbered, insert or update numbered comet record
                # Insert new comet entries from newly downloaded JPL catalog
                mcd.output_log_entry(path_logfile,'Adding/updating data for comet {:s}/{:s}...'.format(comet_number,comet_desig))
                query = "INSERT OR IGNORE INTO ssobjects(object_type,desig_number,desig_name,epoch_mjd,perihelion_dist,eccentricity,inclination,arg_perihelion,long_ascnode,t_perihelion,date_added) VALUES ('comet','{:s}','{:<40s}',{:f},{:f},{:f},{:f},{:f},{:f},{:f},'{:s}')".format(comet_number,
                    comet_desig,epoch,peridist,eccen,inclin,argperi,longascnode,time_peri,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                # Update existing comet entries from newly downloaded JPL catalog
                query = "UPDATE ssobjects SET object_type='comet',desig_name='{:<40s}',epoch_mjd={:f},perihelion_dist={:f},eccentricity={:f},inclination={:f},arg_perihelion={:f},long_ascnode={:f},t_perihelion={:f},date_updated='{:s}' WHERE desig_number='{:s}'".format(comet_desig,
                    epoch,peridist,eccen,inclin,argperi,longascnode,time_peri,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),comet_number)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
            else:  # If input comet is not numbered, insert or update provisional comet record
                # Check if comet entry exists with same provisional designation
                mcd.output_log_entry(path_logfile,'Adding/updating data for comet {:s}...'.format(comet_desig))
                query = "SELECT ssobject_id FROM ssobjects WHERE desig_provisional='{:<40s}' and object_type='comet'".format(comet_desig)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                row = cursor.fetchone()
                if row == None:
                    query = "INSERT OR IGNORE INTO ssobjects(object_type,desig_provisional,epoch_mjd,perihelion_dist,eccentricity,inclination,arg_perihelion,long_ascnode,t_perihelion,date_added) VALUES ('comet','{:<40s}',{:f},{:f},{:f},{:f},{:f},{:f},{:f},'{:s}')".format(comet_desig,
                        epoch,peridist,eccen,inclin,argperi,longascnode,time_peri,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                else:
                    query = "UPDATE ssobjects SET object_type='comet',epoch_mjd={:f},perihelion_dist={:f},eccentricity={:f},inclination={:f},arg_perihelion={:f},long_ascnode={:f},t_perihelion={:f},date_updated='{:s}' WHERE desig_provisional='{:<40s}'".format(epoch,
                        peridist,eccen,inclin,argperi,longascnode,time_peri,datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'),comet_desig)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
    # Committing changes and closing the connection to the database file
    conn.commit()
    conn.close()
    mcd.output_log_entry(path_logfile,'Insertion of osculating orbital elements for comets from {:s} complete.'.format(comet_filepath))
    return None


def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_add_update_orbelements.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DB03_add_update_orbelements execution started','DB03_add_update_orbelements execution started.')
    
        current_date   = datetime.datetime.today().strftime('%Y%m%d')
        dir_downloads  = base_path + 'orbital_elements/'

        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DB03_add_update_orbelements')
                
        if not os.path.isdir(dir_downloads): os.mkdir(dir_downloads)

        numbr_url      = 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR'
        unnum_url      = 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM'
        comet_url      = 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET'
        numbr_filename = 'ELEMENTS_NUMBR_{:s}.dat'.format(current_date)
        unnum_filename = 'ELEMENTS_UNNUM_{:s}.dat'.format(current_date)
        comet_filename = 'ELEMENTS_COMET_{:s}.dat'.format(current_date)
            
        numbr_filepath = dir_downloads + numbr_filename
        unnum_filepath = dir_downloads + unnum_filename
        comet_filepath = dir_downloads + comet_filename

        # Connect to database file
        conn = mcd.create_connection(sqlite_file)
        cursor = conn.cursor()

        download_jpl_orbelem_files(numbr_url,unnum_url,comet_url,numbr_filepath,unnum_filepath,comet_filepath,path_logfile,path_errorfile)

        ##### For testing only #####
        #numbr_filename = 'ELEMENTS_NUMBR_testing.dat'.format(current_date)
        #numbr_filepath = dir_downloads + numbr_filename

        ingest_numbr_data(numbr_filepath,sqlite_file,path_logfile,path_errorfile)
        #ingest_unnum_data(unnum_filepath,sqlite_file,path_logfile,path_errorfile)
        #ingest_comet_data(comet_filepath,sqlite_file,path_logfile,path_errorfile)

        # Committing changes and closing the connection to the database file
        conn.commit()
        conn.close()
        mcd.output_log_entry(path_logfile,'Database closed.')
    
        mcd.send_status_email('DB03_add_update_orbelements execution complete','DB03_add_update_orbelements execution started.')
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None


if __name__ == '__main__':
    main()

