import re
import sys
import time
import datetime
import math
import numpy as np
from numpy import linspace
import glob, os, bz2, subprocess
import os.path
import sqlite3
from sqlite3 import Error
import datetime
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- SEARCH HISTORY INGESTION FUNCTIONS #####

def ingest_search_history_file(search_history_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile):
    try:
        with open(search_history_file_to_ingest) as input_file:
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            for _ in range(1): #skip one header line
                next(input_file)
            for line in input_file:
                ssobject_id        = int(line.split()[0])
                instrument_id      = int(line.split()[1])
                last_searched_date = line.split()[2]
                query = "SELECT ssobject_id,instrument_id,last_searched FROM search_history WHERE ssobject_id={:d} AND instrument_id={:d}".format(ssobject_id,instrument_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                row_detection = cursor.fetchone()
                if row_detection == None:
                    query = "INSERT INTO search_history(ssobject_id,instrument_id,last_searched) VALUES ({:d},{:d},'{:s}')".format(ssobject_id,instrument_id,last_searched_date)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                else:
                    query = "UPDATE search_history SET last_searched='{:s}' WHERE ssobject_id={:d} AND instrument_id={:d}".format(last_searched_date,ssobject_id,instrument_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
            conn.commit() # Commit changes
            conn.close()  # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: ingest_search_history_file()')
        keep_going = False
        print(e)
    return keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_ingest_exposure_data_sdss.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('SQ02_ingest_search_history execution started','{:s} - SQ02_ingest_search_history execution started'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'SQ02_ingest_search_history')

        ssois_directory = base_path + 'ssois_queries/'
        search_urls_directory = ssois_directory + 'search_urls/'
        if not os.path.isdir(ssois_directory):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS query directory {:s} not found'.format(ssois_directory))
            mcd.send_status_email('SQ02_ingest_search_history execution failed','SQ02_ingest_search_history execution failed - SSOIS query directory {:s} not found'.format(ssois_directory))
            keep_going = False
        if not os.path.isdir(search_urls_directory):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS search URLs directory {:s} not found'.format(search_urls_directory))
            mcd.send_status_email('SQ02_ingest_search_history execution failed','{:s} - SQ02_ingest_search_history execution failed - SSOIS search URLs directory {:s} not found'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),search_urls_directory))
            keep_going = False

    if keep_going:
        # Connect to database file
        os.chdir(search_urls_directory)
        mcd.output_log_entry(path_logfile,'Ingesting SSOIS search results into {:s}...'.format(sqlite_file))
        for search_history_file_to_ingest_gz in sorted(glob.glob('search_history_*_toingest.txt.gz')):
            mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(search_history_file_to_ingest_gz))
            mcd.decompress_file_gzip(search_history_file_to_ingest_gz)
            search_history_file_to_ingest = search_history_file_to_ingest_gz[:-3]
            ingest_search_history_file(search_history_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile)
            if keep_going:
                ingested_filename = search_history_file_to_ingest[:-13]+'_ingested.txt'
                os.rename(search_history_file_to_ingest,ingested_filename)
                mcd.compress_file_gzip(ingested_filename)
            else:
                mcd.compress_file_gzip(search_history_file_to_ingest)

    mcd.send_status_email('SQ02_ingest_search_history execution complete','SQ02_ingest_search_history execution complete')
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()
