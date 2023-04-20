import sys
import glob, os, bz2, subprocess
import os.path
import datetime
import sqlite3
from sqlite3 import Error
import smtplib, ssl
import macadamia_functions as mcd

### Edit History ###
# 2021-06-16: implemented new logging functionality (new logging functions and more frequent log updates)
# Updated 2021-06-16, 20:50

##### FUNCTION DEFINITIONS -- DATABASE OPERATIONS #####

def ingest_photometry_data(sqlite_file,keep_going,path_logfile,path_errorfile):
    if keep_going:
        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
        cursor = conn.cursor()
        for photometry_data_to_ingest_filename_gz in sorted(glob.glob('00_photometry_data_*_toingest.txt.gz')):
            try:
                mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(photometry_data_to_ingest_filename_gz))
                mcd.decompress_file_gzip(photometry_data_to_ingest_filename_gz)
                photometry_data_to_ingest_filename = photometry_data_to_ingest_filename_gz[:-3]                
                with open(photometry_data_to_ingest_filename,'r') as expdata_toingest:
                    for _ in range(1): # Skip first header line
                        next(expdata_toingest)
                    for line in expdata_toingest:
                        expdata = line.split()
                        mosaic_element_id     = int(expdata[0])
                        base_filename         = expdata[1]
                        date_tai              = expdata[2]
                        exposure_start_tai    = expdata[3]
                        filter_name           = expdata[4]
                        if expdata[5] != 'nan': pointing_center_ra  = float(expdata[5])
                        else:                   pointing_center_ra  = -99
                        if expdata[6] != 'nan': pointing_center_dec = float(expdata[6])
                        else:                   pointing_center_dec = -99
                        sky_mean              = float(expdata[7])
                        sky_stddev            = float(expdata[8])
                        psf_width_mean_arcsec = float(expdata[9])
                        source_density        = float(expdata[10])
                        avg_zeropoint_mag     = float(expdata[11])
                        avg_zeropoint_err     = float(expdata[12])
                        num_calib_stars       = int(expdata[13])
                        limit_mag_ps          = float(expdata[14])
                        limit_mag_sb          = float(expdata[15])
                        source_table_path     = expdata[16]
                        source_table_file     = expdata[17]

                        # Check for existing exposure ID with matching mosaic element, date, time, filter, and pointing center
                        query = "SELECT exposure_id FROM exposures WHERE mosaic_element_id={:d} AND base_filename='{:s}'".format(mosaic_element_id,base_filename)
                        mcd.output_log_entry(path_logfile,query)
                        cursor.execute(query)
                        row = cursor.fetchone()
                        if row == None: # If exposure entry not found
                            mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure database entry for mosaic_element_id={:d} and base_filename={:s} not found.'.format(mosaic_element_id,base_filename))
                        else:
                            mcd.output_log_entry(path_logfile,'Updating exposure database entry.')
                            exposure_id = int(row[0])
                            query = "UPDATE exposures SET source_table_path='{:s}',source_table_file='{:s}',sky_mean={:.6f},sky_stddev={:.6f},psf_width_mean={:.3f},source_density={:.3f},zero_point={:.3f},zpoint_err={:.3f},zpoint_nstars={:d},limiting_mag_ps_3sig={:.3f},limiting_mag_sb_3sig={:.3f} WHERE exposure_id={:d}".format(source_table_path,source_table_file,sky_mean,sky_stddev,psf_width_mean_arcsec,source_density,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,limit_mag_ps,limit_mag_sb,exposure_id)
                            mcd.output_log_entry(path_logfile,query)
                            cursor.execute(query)
                photometry_data_ingested_filename = photometry_data_to_ingest_filename[:-13]+'_ingested.txt'
                mcd.rename_file(photometry_data_to_ingest_filename,photometry_data_ingested_filename)
                mcd.compress_file_gzip(photometry_data_ingested_filename)
            except Exception as e:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: ingest_photometry_data_sdss()'.format(photometry_data_to_ingest_filename_gz))
                print(e)
        conn.commit() # Commit changes
        conn.close()  # Close connection to database file
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: ingest_photometry_data_sdss()')
    return keep_going


def insert_exp_photometry_status(sqlite_file,keep_going,path_logfile,path_errorfile):
    if keep_going:
        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
        cursor = conn.cursor()
        for processed_files_toingest_filename_gz in sorted(glob.glob('00_processed_files_photometry_*_toingest.txt.gz')):
            try:
                mcd.decompress_file_gzip(processed_files_toingest_filename_gz)
                processed_files_toingest_filename = processed_files_toingest_filename_gz[:-3]
                with open(processed_files_toingest_filename,'r') as inputfile:
                    for _ in range(1): # Skip first header line
                        next(inputfile)
                    for line in inputfile:
                        fits_filename_fz    = line[:50].strip()
                        status              = line[53:-1].strip()
                        proc_data_file      = fits_filename_fz
                        query = "UPDATE exposures SET exposure_status='{:s}' WHERE proc_data_file='{:s}'".format(status,proc_data_file)
                        mcd.output_log_entry(path_logfile,query)
                        cursor.execute(query)
                processed_files_ingested_filename = processed_files_toingest_filename[:-13] + '_ingested.txt'
                mcd.rename_file(processed_files_toingest_filename,processed_files_ingested_filename)
                mcd.compress_file_gzip(processed_files_ingested_filename)
            except Exception as e:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: insert_exp_status()'.format(processed_files_toingest_filename_gz))
                print(e)
        conn.commit() # Commit changes
        conn.close()  # Close connection to database file
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: insert_exp_status()')
    return keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_IM05b_ingest_photometry_data.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('IM05b_ingest_photometry_data execution started','IM05b_ingest_photometry_data execution started.')

        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'IM05b_ingest_photometry_data')
        
        os.chdir(base_path)
        dir_exposure_data = 'exposure_data/'
        if not os.path.isdir(dir_exposure_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure data directory {:s} not found.'.format(dir_exposure_data))
            mcd.send_status_email('IM05b_ingest_photometry_data execution failed','IM05b_ingest_photometry_data execution failed - Exposure data directory {:s} not found.'.format(dir_exposure_data))
            keep_going = False

    if keep_going:
        # Ingest data
        os.chdir(base_path+dir_exposure_data)
        keep_going = ingest_photometry_data(sqlite_file,keep_going,path_logfile,path_errorfile)
        keep_going = insert_exp_photometry_status(sqlite_file,keep_going,path_logfile,path_errorfile)
            
    mcd.send_status_email('IM05b_ingest_photometry_data execution complete','IM05b_ingest_photometry_data execution complete.')
                
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None


if __name__ == '__main__':
    main()

