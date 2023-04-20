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

def ingest_detection_details_file(detection_details_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile):
    try:
        with open(detection_details_file_to_ingest) as input_file:
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            for _ in range(1): #skip one header line
                next(input_file)
            for line in input_file:
                det_details             = line.split()
                detection_id            = int(det_details[0])
                ra_predicted            = float(det_details[1])
                dec_predicted           = float(det_details[2])
                ra_predicted_err        = float(det_details[3])
                dec_predicted_err       = float(det_details[4])
                ra_rate                 = float(det_details[5])
                dec_rate                = float(det_details[6])
                trail_pa_pred           = float(det_details[7])
                trail_len_pred          = float(det_details[8])
                phase_angle             = float(det_details[9])
                solar_elongation        = float(det_details[10])
                heliocentric_dist       = float(det_details[11])
                geocentric_dist         = float(det_details[12])
                vmag_predicted          = float(det_details[13])
                pa_antisolar            = float(det_details[14])
                pa_neg_heliocentric_vel = float(det_details[15])
                orbit_plane_angle       = float(det_details[16])
                ecliptic_longitude      = float(det_details[17])
                ecliptic_latitude       = float(det_details[18])
                lunar_illumination      = float(det_details[19])
                lunar_elongation        = float(det_details[20])
                true_anomaly            = float(det_details[21])
                galactic_longitude      = float(det_details[22])
                galactic_latitude       = float(det_details[23])
            
                query = "SELECT detection_id FROM detection_details WHERE detection_id={:d}".format(detection_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                row_detection = cursor.fetchone()
                if row_detection == None:
                    mcd.output_log_entry(path_logfile,'Inserting detection_details entry for detection_id {:d}...'.format(detection_id))
                    query = "INSERT OR IGNORE INTO detection_details(detection_id,ra_predicted,dec_predicted,ra_predicted_err,dec_predicted_err,ra_rate,dec_rate,trail_pa_expected,trail_length_expected,phase_angle,solar_elongation,heliocentric_dist,geocentric_dist,vmag_predicted,pa_antisolar,pa_neg_heliocentric_vel,orb_plane_angle,ecliptic_longitude,ecliptic_latitude,lunar_illumination,lunar_elongation,true_anomaly,galactic_longitude,galactic_latitude,detection_details_status) VALUES ({:d},{:.7f},{:.7f},{:.7f},{:.7f},{:.7f},{:.7f},{:.2f},{:.2f},{:.4f},{:.4f},{:.8f},{:.8f},{:.2f},{:.3f},{:.3f},{:.4f},{:.4f},{:.4f},{:.2f},{:.2f},{:.4f},{:.4f},{:.4f},'Added, {:s}')".format(detection_id,ra_predicted,dec_predicted,ra_predicted_err,dec_predicted_err,ra_rate,dec_rate,trail_pa_pred,trail_len_pred,phase_angle,solar_elongation,heliocentric_dist,geocentric_dist,vmag_predicted,pa_antisolar,pa_neg_heliocentric_vel,orbit_plane_angle,ecliptic_longitude,ecliptic_latitude,lunar_illumination,lunar_elongation,true_anomaly,galactic_longitude,galactic_latitude,datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                else:
                    mcd.output_log_entry(path_logfile,'Updating detection_details entry for detection_id {:d}...'.format(detection_id))
                    query = "UPDATE detection_details SET ra_predicted={:.7f},dec_predicted={:.7f},ra_predicted_err={:.7f},dec_predicted_err={:.7f},ra_rate={:.7f},dec_rate={:.7f},trail_pa_expected={:.2f},trail_length_expected={:.2f},phase_angle={:.4f},solar_elongation={:.4f},heliocentric_dist={:.8f},geocentric_dist={:.8f},vmag_predicted={:.2f},pa_antisolar={:.3f},pa_neg_heliocentric_vel={:.3f},orb_plane_angle={:.3f},ecliptic_longitude={:.4f},ecliptic_latitude={:.4f},lunar_illumination={:.2f},lunar_elongation={:.2f},true_anomaly={:.4f},galactic_longitude={:.4f},galactic_latitude={:.4f},detection_details_status='Updated, {:s}' WHERE detection_id={:d}".format(ra_predicted,dec_predicted,ra_predicted_err,dec_predicted_err,ra_rate,dec_rate,trail_pa_pred,trail_len_pred,phase_angle,solar_elongation,heliocentric_dist,geocentric_dist,vmag_predicted,pa_antisolar,pa_neg_heliocentric_vel,orbit_plane_angle,ecliptic_longitude,ecliptic_latitude,lunar_illumination,lunar_elongation,true_anomaly,galactic_longitude,galactic_latitude,datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),detection_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
            conn.commit() # Commit changes
            conn.close()  # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: ingest_detection_details_file'.format(detection_details_file_to_ingest))
        keep_going = False
        print(e)
    return keep_going


def main():
    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_ingest_detection_exposure_object_links.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DT05_ingest_detection_details execution started','DT05_ingest_detection_details execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DT05_ingest_detection_details')

        dir_detection_data = base_path + 'detection_data/'
        if not os.path.isdir(dir_detection_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Detection data directory {:s} not found'.format(dir_detection_data))
            mcd.send_status_email('DT05_ingest_detection_details execution failed','DT05_ingest_detection_details execution failed - Detection data directory {:s} not found.'.format(dir_detection_data))
            keep_going = False

    if keep_going:
        os.chdir(dir_detection_data)
        for detection_details_file_to_ingest_gz in sorted(glob.glob('detection_details_*_toingest.txt.gz')):
            mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(detection_details_file_to_ingest_gz))
            mcd.decompress_file_gzip(detection_details_file_to_ingest_gz)
            detection_details_file_to_ingest = detection_details_file_to_ingest_gz[:-3]
            keep_going = ingest_detection_details_file(detection_details_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile)
            if keep_going:
                ingested_filename = detection_details_file_to_ingest[:-13]+'_ingested.txt'
                os.rename(detection_details_file_to_ingest,ingested_filename)
                mcd.compress_file_gzip(ingested_filename)
            else:
                mcd.compress_file_gzip(detection_details_file_to_ingest)

    mcd.send_status_email('DT05_ingest_detection_details execution complete','DT05_ingest_detection_details execution complete.')
    
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

