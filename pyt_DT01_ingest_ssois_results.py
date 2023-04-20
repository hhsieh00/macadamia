import re
import sys
import time
import datetime
import math
import urllib
import glob, os, bz2, subprocess
import os.path
import sqlite3
from sqlite3 import Error
import requests
import warnings
import statistics
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE FUNCTIONS #####

def ingest_results_file(results_file,telinst,sqlite_file,path_logfile,path_errorfile):
    conn = mcd.create_connection(sqlite_file)
    cursor = conn.cursor()
    with open(results_file) as input_file:
        for _ in range(1): #skip one header line
            next(input_file)
        for line in input_file:
            search_result_object  = mcd.add_apostrophe_escape(line[0:10].strip())
            search_result_date    = line[12:22].strip()
            search_result_time    = line[23:31].strip()
            search_result_filter  = line[37:47].strip()
            search_result_exptime = float(line[49:56])
            search_result_ra      = float(line[58:68])
            search_result_dec     = float(line[70:80])
            search_result_target  = line[82:112].strip()
            search_result_telinst = line[114:134].strip()
            search_result_datalink = line[136:].strip()
            if search_result_filter[:4] == 'G.MP': search_result_filter = 'g'
            if search_result_filter[:4] == 'R.MP': search_result_filter = 'r'
            if search_result_filter[:4] == 'I.MP': search_result_filter = 'i'
            if search_result_filter[:4] == 'Z.MP': search_result_filter = 'z'            
            if search_result_telinst == telinst and search_result_filter != 'u':
                mcd.output_log_entry(path_logfile,'Adding search result entry to database ({:s},{:s},{:s},{:s},{:s})...'.format(search_result_object,search_result_date,search_result_time,search_result_telinst,search_result_filter))
                query = "INSERT OR IGNORE INTO detections(search_result_object,search_result_date,search_result_time,search_result_filter,search_result_exptime,search_result_ra,search_result_dec,search_result_target,search_result_telinst,search_result_datalink,search_result_status) VALUES ('{:s}','{:s}','{:s}','{:s}',{:.3f},{:.6f},{:.6f},'{:s}','{:s}','{:s}','Detection ingested.')".format(search_result_object,search_result_date,search_result_time,search_result_filter,search_result_exptime,search_result_ra,search_result_dec,search_result_target,search_result_telinst,search_result_datalink)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
    conn.commit() # Commit changes
    conn.close()  # Close connection to database file
    return None


def main():

    # Define filenames and paths
    if len(sys.argv)!=4:
        print('Usage:\n python3 pyt_DT01_ingest_ssois_results.py [base_path] [sqlite_file] [Telescope/Instrument]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    telinst     = sys.argv[3]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DT01_ingest_ssois_results execution started','DT01_ingest_ssois_results execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DT01_ingest_ssois_results')

        ssois_directory = base_path + 'ssois_queries/parsed_results/'
        if not os.path.isdir(ssois_directory):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Parsed SSOIS results directory {:s} not found'.format(ssois_directory))
            mcd.send_status_email('DT01_ingest_ssois_results execution failed','DT01_ingest_ssois_results execution failed - Parsed SSOIS results directory {:s} not found.'.format(ssois_directory))
            keep_going = False

    if keep_going:
        os.chdir(ssois_directory)
        for results_file_to_ingest_gz in sorted(glob.glob('parsed_ssois_results_*_toingest.txt.gz')):
            mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(results_file_to_ingest_gz))
            mcd.decompress_file_gzip(results_file_to_ingest_gz)
            results_file_to_ingest = results_file_to_ingest_gz[:-3]
            ingest_results_file(results_file_to_ingest,telinst,sqlite_file,path_logfile,path_errorfile)
            ingested_filename = results_file_to_ingest[:-13]+'_ingested.txt'
            os.rename(results_file_to_ingest,ingested_filename)
            mcd.compress_file_gzip(ingested_filename)

    mcd.send_status_email('DT01_ingest_ssois_results execution complete','DT01_ingest_ssois_results execution complete.')
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

