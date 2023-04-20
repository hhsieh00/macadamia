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
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE OPERATIONS #####

def ingest_det_exp_obj_file(det_exp_obj_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile):
    try:
        with open(det_exp_obj_file_to_ingest) as input_file:
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            for _ in range(1): #skip one header line
                next(input_file)
            for line in input_file:
                det_exp_obj  = line.split()
                detection_id = int(det_exp_obj[0])
                exposure_id  = int(det_exp_obj[1])
                ssobject_id  = int(det_exp_obj[2])
                if exposure_id == -1:
                    status = 'Detection not in FOV of specified exposure.'
                else:
                    status = 'Exposure, object, and detection linked.'

                query = "SELECT detection_id FROM detection_data WHERE detection_id={:d}".format(detection_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                row_detection = cursor.fetchone()
                if row_detection == None:
                    mcd.output_log_entry(path_logfile,'Inserting detection_data entry for detection_id {:d}...'.format(detection_id))
                    query = "INSERT OR IGNORE INTO detection_data(detection_id,exposure_id,ssobject_id,detection_status) VALUES ({:d},{:d},{:d},'{:s}')".format(detection_id,exposure_id,ssobject_id,status)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                    query = "UPDATE detections SET search_result_status='{:s}' WHERE detection_id={:d}".format(status,detection_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                else:
                    mcd.output_log_entry(path_logfile,'Updating detections entry for detection_id {:d}...'.format(detection_id))
                    query = "UPDATE detection_data SET exposure_id={:d},ssobject_id={:d},detection_status='{:s}' WHERE detection_id={:d}".format(exposure_id,ssobject_id,status,detection_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                    query = "UPDATE detections SET search_result_status='{:s}' WHERE detection_id={:d}".format(status,detection_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
            conn.commit() # Commit changes
            conn.close()  # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: ingest_det_exp_obj_file()'.format(det_exp_obj_file_to_ingest))
        keep_going = False
        print(e)
    return keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_DT03_ingest_detection_exposure_object_links.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DT03_ingest_detection_exposure_object_links execution started','DT03_ingest_detection_exposure_object_links execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DT03_ingest_detection_exposure_object_links')

        dir_detection_data = base_path + 'detection_data/'
        if not os.path.isdir(dir_detection_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Detection data directory {:s} not found'.format(dir_detection_data))
            mcd.send_status_email('DT03_ingest_detection_exposure_object_links execution failed','DT03_ingest_detection_exposure_object_links execution failed - Detection data directory {:s} not found.'.format(dir_detection_data))
            keep_going = False

    if keep_going:
        os.chdir(dir_detection_data)
        for det_exp_obj_file_to_ingest_gz in sorted(glob.glob('detection_exposure_object_links_*_toingest.txt.gz')):
            mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(det_exp_obj_file_to_ingest_gz))
            mcd.decompress_file_gzip(det_exp_obj_file_to_ingest_gz)
            det_exp_obj_file_to_ingest = det_exp_obj_file_to_ingest_gz[:-3]
            keep_going = ingest_det_exp_obj_file(det_exp_obj_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile)
            if keep_going:
                ingested_filename = det_exp_obj_file_to_ingest[:-13]+'_ingested.txt'
                os.rename(det_exp_obj_file_to_ingest,ingested_filename)
                mcd.output_log_entry(path_logfile,'Compressing {:s}...'.format(ingested_filename))
                mcd.compress_file_gzip(ingested_filename)
            else:
                mcd.output_log_entry(path_logfile,'Compressing {:s}...'.format(det_exp_obj_file_to_ingest))
                mcd.compress_file_gzip(det_exp_obj_file_to_ingest)
                
    mcd.send_status_email('DT03_ingest_detection_exposure_object_links execution complete','{:s} - DT03_ingest_detection_exposure_object_links execution complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

