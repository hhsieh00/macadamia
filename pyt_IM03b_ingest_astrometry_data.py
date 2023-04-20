import sys
import glob, os, bz2, subprocess
import os.path
import datetime
import sqlite3
from sqlite3 import Error
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE OPERATIONS #####

def ingest_astrometry_data(sqlite_file,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            for astrometry_data_to_ingest_filename_gz in sorted(glob.glob('00_astrometry_data_*_toingest.txt.gz')):
                mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(astrometry_data_to_ingest_filename_gz))
                
                mcd.decompress_file_gzip(astrometry_data_to_ingest_filename_gz)
                astrometry_data_to_ingest_filename = astrometry_data_to_ingest_filename_gz[:-3]                
                with open(astrometry_data_to_ingest_filename) as expdata_toingest:
                    for _ in range(1): #skip first header line
                        next(expdata_toingest)
                    for line in expdata_toingest:
                        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
                        cursor = conn.cursor()
                
                        expdata = line.split()
                        mosaic_element_id     = int(expdata[0])
                        base_filename         = expdata[1]
                        date_tai              = expdata[2]
                        exposure_start_jd     = float(expdata[3])
                        exposure_start_tai    = expdata[4]
                        exposure_time         = float(expdata[5])
                        filter_name           = expdata[6]
                        
                        airmass_str             = expdata[7]
                        if airmass_str == 'nan' or airmass_str == 'inf':
                            airmass_str = '0.000'
                        airmass                 = float(airmass_str)
                        
                        pointing_center_ra_str  = expdata[8]
                        if pointing_center_ra_str == 'nan' or pointing_center_ra_str == 'inf':
                            pointing_center_ra_str = '0.000'
                        pointing_center_ra      = float(pointing_center_ra_str)
                        
                        pointing_center_dec_str = expdata[9]
                        if pointing_center_dec_str == 'nan' or pointing_center_dec_str == 'inf':
                            pointing_center_dec_str = '0.000'
                        pointing_center_dec     = float(pointing_center_dec_str)
                        
                        tracking_rate_ra      = float(expdata[10])
                        tracking_rate_dec     = float(expdata[11])
                        binning_x             = int(expdata[12])
                        binning_y             = int(expdata[13])
                        exp_npix_x            = int(expdata[14])
                        exp_npix_y            = int(expdata[15])
                        wcs_nstars            = int(expdata[16])
                        wcs_crpix1            = int(expdata[17])
                        wcs_crpix2            = int(expdata[18])
                        wcs_crval1            = float(expdata[19])
                        wcs_crval2            = float(expdata[20])
                        wcs_cd1_1             = float(expdata[21])
                        wcs_cd1_2             = float(expdata[22])
                        wcs_cd2_1             = float(expdata[23])
                        wcs_cd2_2             = float(expdata[24])
                        wcs_pixscale_x        = float(expdata[25])
                        wcs_pixscale_y        = float(expdata[26])
                        raw_data_link         = expdata[27]
                        raw_data_path         = expdata[28]
                        raw_data_file         = expdata[29]
                        proc_data_path        = expdata[30]
                        proc_data_file        = expdata[31]

                        # Check for existing exposure ID with matching mosaic element, date, time, filter, and pointing center
                        #query = "SELECT exposure_id FROM exposures WHERE mosaic_element_id={:d} AND date_tai='{:s}' AND exposure_start_tai='{:s}' AND filter_name='{:s}' AND pointing_center_ra={:f} AND pointing_center_dec={:f} AND base_filename='{:s}'".format(mosaic_element_id,date_tai,exposure_start_tai,filter_name,pointing_center_ra,pointing_center_dec,base_filename)
                        query = "SELECT exposure_id FROM exposures WHERE mosaic_element_id={:d} AND base_filename='{:s}'".format(mosaic_element_id,base_filename)
                        mcd.output_log_entry(path_logfile,query)
                        cursor.execute(query)
                        row = cursor.fetchone()
                        if row == None: # If no exposure ID entry exists, add new entry
                            mcd.output_log_entry(path_logfile,'Creating new exposure database entry.')
                            query = "INSERT INTO exposures(mosaic_element_id,base_filename,raw_data_link,raw_data_path,raw_data_file,proc_data_path,proc_data_file,date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra,pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exposure_npix_x,exposure_npix_y,wcs_nstars,wcs_crpix1,wcs_crpix2,wcs_crval1,wcs_crval2,wcs_cd1_1,wcs_cd1_2,wcs_cd2_1,wcs_cd2_2,wcs_pixscale_x,wcs_pixscale_y) VALUES ({:d},'{:s}','{:s}','{:s}','{:s}','{:s}','{:s}','{:s}',{:.6f},'{:s}',{:.6f},'{:s}',{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:d},{:d},{:d},{:d},{:d},{:d},{:d},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f})".format(mosaic_element_id,base_filename,raw_data_link,raw_data_path,raw_data_file,proc_data_path,proc_data_file,date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra,pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,wcs_nstars,wcs_crpix1,wcs_crpix2,wcs_crval1,wcs_crval2,wcs_cd1_1,wcs_cd1_2,wcs_cd2_1,wcs_cd2_2,wcs_pixscale_x,wcs_pixscale_y)
                            mcd.output_log_entry(path_logfile,query)
                            cursor.execute(query)
                            exposure_id = int(cursor.lastrowid)
                        else:
                            exposure_id = int(row[0])
                            query = "UPDATE exposures SET raw_data_link='{:s}',raw_data_path='{:s}',raw_data_file='{:s}',proc_data_path='{:s}',proc_data_file='{:s}',exposure_start_jd={:.10f},exposure_time={:.6f},airmass={:.6f},tracking_rate_ra={:.6f},tracking_rate_dec={:.6f},binning_x={:d},binning_y={:d},exposure_npix_x={:d},exposure_npix_y={:d},wcs_nstars={:d},wcs_crpix1={:d},wcs_crpix2={:d},wcs_crval1={:.6f},wcs_crval2={:.6f},wcs_cd1_1={:.6f},wcs_cd1_2={:.6f},wcs_cd2_1={:.6f},wcs_cd2_2={:.6f},wcs_pixscale_x={:.6f},wcs_pixscale_y={:.6f} WHERE exposure_id={:d}".format(raw_data_link,raw_data_path,raw_data_file,proc_data_path,proc_data_file,exposure_start_jd,exposure_time,airmass,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,wcs_nstars,wcs_crpix1,wcs_crpix2,wcs_crval1,wcs_crval2,wcs_cd1_1,wcs_cd1_2,wcs_cd2_1,wcs_cd2_2,wcs_pixscale_x,wcs_pixscale_y,exposure_id)
                            mcd.output_log_entry(path_logfile,query)
                            cursor.execute(query)
                            mcd.output_log_entry(path_logfile,'Updating exposure database entry.')
                
                        conn.commit() # Commit changes
                        conn.close()  # Close connection to database file
                        
                astrometry_data_ingested_filename = astrometry_data_to_ingest_filename[:-13] + '_ingested.txt'
                os.rename(astrometry_data_to_ingest_filename,astrometry_data_ingested_filename)
                mcd.compress_file_gzip(astrometry_data_ingested_filename)
                
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: ingest_astrometry_data()')
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: ingest_astrometry_data()')
    return keep_going


def insert_exp_astrometry_status(sqlite_file,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            for processed_files_toingest_filename_gz in sorted(glob.glob('00_processed_files_astrometry_*_toingest.txt.gz')):
                mcd.decompress_file_gzip(processed_files_toingest_filename_gz)
                processed_files_toingest_filename = processed_files_toingest_filename_gz[:-3]
                with open(processed_files_toingest_filename,'r') as inputfile:
                    for _ in range(1): #skip first header line
                        next(inputfile)
                    for line in inputfile:
                        fits_filename_fz    = line[:50].strip()
                        status              = line[53:-1].strip()
                        #fits_filestem       = fits_filename_fz[:-12]  # e.g., frame-z-004828-1-0430_000_wcs.fits.fz --> frame-z-004828-1-0430_000
                        proc_data_file      = fits_filename_fz[:-8] + '_wcs.fits.fz'  # e.g., frame-z-004828-1-0430_000.fits.fz --> frame-z-004828-1-0430_000_wcs.fits.fz
                        update_status_query = "UPDATE exposures SET exposure_status='{:s}' WHERE proc_data_file='{:s}'".format(status,proc_data_file)
                        mcd.output_log_entry(path_logfile,update_status_query)
                        cursor.execute(update_status_query)
                processed_files_ingested_filename = processed_files_toingest_filename[:-13] + '_ingested.txt'
                mcd.rename_file(processed_files_toingest_filename,processed_files_ingested_filename)
                mcd.compress_file_gzip(processed_files_ingested_filename)
            conn.commit() # Commit changes
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: insert_exp_status()'.format(dir_processed_data))
            mcd.output_error_log_entry(path_logfile,path_errorfile,update_status_query)
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: insert_exp_status()')
    return keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_IM03b_ingest_exposure_data.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('IM03b_ingest_astrometry_data execution started','IM03b_ingest_astrometry_data execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'IM03b_ingest_astrometry_data')
        
        #os.chdir(base_path)
        dir_processed_data = 'exposure_data/'
        if not os.path.isdir(base_path + dir_processed_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Processed data directory {:s} not found.'.format(dir_processed_data))
            mcd.send_status_email('IM03b_ingest_astrometry_data execution failed','IM03b_ingest_astrometry_data execution failed - Processed data directory {:s} not found.'.format(dir_processed_data))
            keep_going = False

    if keep_going:
        # Ingest data
        os.chdir(base_path+dir_processed_data)
        keep_going = ingest_astrometry_data(sqlite_file,keep_going,path_logfile,path_errorfile)
        keep_going = insert_exp_astrometry_status(sqlite_file,keep_going,path_logfile,path_errorfile)

        mcd.send_status_email('IM03b_ingest_astrometry_data execution complete','IM03b_ingest_astrometry_data execution complete.')
                
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
            
    return None


if __name__ == '__main__':
    main()

