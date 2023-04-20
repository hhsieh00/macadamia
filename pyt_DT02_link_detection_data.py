import re
import sys
import datetime
import math
import numpy as np
import glob, os, bz2, subprocess
import os.path
import sqlite3
from sqlite3 import Error
from astroquery.gaia import Gaia
from astroquery.jplhorizons import Horizons
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DETECTION-OBJECT-EXPOSURE LINKING #####

def get_filename_from_datalink(datalink):
    num_chars = len(datalink)
    idx = num_chars - 1
    while datalink[idx] != '/':
        idx -= 1
    filename = datalink[(idx+1):]
    return filename


def link_detections_exposures_ssobjects(sqlite_file,output_file_links,exposures_not_found_file,instrument,thread_idx,start_detection_id,path_logfile,path_errorfile):
    try:
        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
        cursor = conn.cursor()

        mcd.output_log_entry(path_logfile,'Retrieving unlinked detections...')
        if instrument == 'SDSS':
            query = "SELECT detection_id,search_result_object,search_result_date,search_result_telinst,search_result_datalink FROM detections WHERE search_result_status='Detection ingested.' AND search_result_telinst='SDSS'"
        elif instrument == 'MegaCam':
            query = "SELECT detection_id,search_result_object,search_result_date,search_result_telinst,search_result_datalink FROM detections WHERE search_result_status='Detection ingested.' AND search_result_telinst='CFHT/MegaCam'"
        else:
            query = "SELECT detection_id,search_result_object,search_result_date,search_result_telinst,search_result_datalink FROM detections WHERE search_result_status='Detection ingested.'"
        mcd.output_log_entry(path_logfile,query)
        cursor.execute(query)
        rows = cursor.fetchall()
        mcd.output_log_entry(path_logfile,'{:d} unlinked detections retrieved'.format(len(rows)))
        for row in rows:
            keep_going = True
            detection_id           = row[0]
            search_result_object   = row[1]
            search_result_date     = row[2]
            search_result_telinst  = row[3]
            search_result_datalink = row[4]
            
            if '{:010d}'.format(detection_id)[-2:] == '{:02d}'.format(thread_idx) and detection_id > start_detection_id:
                mcd.output_log_entry(path_logfile,'Searching for object and exposure corresponding to Detection {:d}...'.format(detection_id))
            
                ### Retrieve ssobject_id corresponding to given object and retrieve orbital elements
                ssobject_id,desig_number,keep_going = link_detection_ssobject(detection_id,search_result_object,sqlite_file,keep_going,path_logfile,path_errorfile)
            
                ### Retrieve observation time from first available data row corresponding to identified exposure
                ###  and compute expected RA/Dec at time (ignore whether object is in that specific extension for now)
                temp_exposure_id,proc_data_path,temp_data_file,keep_going = link_detection_temp_exposure(detection_id,search_result_date,search_result_telinst,search_result_datalink,exposures_not_found_file,sqlite_file,keep_going,path_logfile,path_errorfile)
                jd_mid,keep_going = compute_jd_mid(temp_exposure_id,sqlite_file,keep_going,path_logfile,path_errorfile)
                right_ascension,declination,keep_going = retrieve_predicted_coords_horizons(desig_number,jd_mid,keep_going,path_logfile,path_errorfile)
                extension_number_with_detection,keep_going = identify_extension(proc_data_path,temp_data_file,right_ascension,declination,sqlite_file,keep_going,path_logfile,path_errorfile)
                exposure_id,keep_going = retrieve_exposure_id(extension_number_with_detection,proc_data_path,temp_data_file,sqlite_file,keep_going,path_logfile,path_errorfile)
                
                if keep_going and ssobject_id != -1:
                    mcd.output_log_entry(path_logfile,'Detection {:d} linked with exposure {:d} and ssobject {:d}'.format(detection_id,exposure_id,ssobject_id))
                    with open(output_file_links,'a') as of:
                        of.write('{:>12d}  {:>11d}  {:>11d}\n'.format(detection_id,exposure_id,ssobject_id))
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Entries in exposures and ssobjects tables for detection_id {:d} not found'.format(detection_id))
        cursor.close() # Close cursor
        conn.close()   # Close connection to database file
    except Exception as e:
        mcd.output_log_entry(path_logfile,'Function failed: link_detections_exposures_ssobjects()')
        print(e)
    return None


def link_detection_ssobject(detection_id,search_result_object,sqlite_file,keep_going,path_logfile,path_errorfile):
    ssobject_id,desig_number = -1,-1
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            if search_result_object.isnumeric():  # For finding ssobject_id's for numbered asteroids
                query = "SELECT ssobject_id,desig_number FROM ssobjects WHERE desig_number='{:>7d}'".format(int(search_result_object))
                cursor1.execute(query)
                row_ssobject = cursor1.fetchone()
                if row_ssobject != None:  # If corresponding solar system object is found, create detection_data entry
                    ssobject_id  = row_ssobject[0]
                    desig_number = row_ssobject[1]
                    mcd.output_log_entry(path_logfile,'Object (ssobject_id={:d}) found for detection_id {:d}'.format(ssobject_id,detection_id))
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Object not found for detection_id {:d}'.format(detection_id))
                    keep_going = False
            else:
                mcd.output_log_entry(path_logfile,'Non-numbered asteroids not being processed at this time')
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for detection {:d}, object {:d}: link_detection_ssobject()'.format(detection_id,search_result_object))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: link_detection_ssobject()')
    return ssobject_id,desig_number,keep_going


def retrieve_orbital_elements(ssobject_id,sqlite_file,keep_going,path_logfile,path_errorfile):
    # Retrieving orbital elements for given solar system object
    orbit_array = np.array([[0,0,0,0,0,0,0,0,0,0,0,0]],dtype=np.double,order='F')
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            mcd.output_log_entry(path_logfile,'Retrieving orbital elements for ssobject {:d} from database...'.format(ssobject_id))
            epoch_mjd,eccentricity,inclination,arg_perihelion,long_ascnode,h_mag,g_param,semimaj_axis,mean_anomaly,perihelion_dist,t_perihelion = -999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999
            query = "SELECT object_type,h_mag,g_param,epoch_mjd,semimaj_axis,perihelion_dist,eccentricity,inclination,arg_perihelion,long_ascnode,mean_anomaly,t_perihelion FROM ssobjects WHERE ssobject_id={:d}".format(ssobject_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                object_type     = row[0]
                epoch_mjd       = float(row[3])
                eccentricity    = float(row[6])
                inclination     = float(row[7])
                arg_perihelion  = float(row[8])
                long_ascnode    = float(row[9])
                if object_type == 'asteroid':
                    h_mag        = float(row[1])
                    g_param      = float(row[2])
                    semimaj_axis = float(row[4])
                    mean_anomaly = float(row[10])
                    if semimaj_axis!=-999 and eccentricity!=-999 and inclination!=-999 and long_ascnode!=-999 and arg_perihelion!=-999 and mean_anomaly!=-999 and epoch_mjd!=-999 and h_mag!=-999 and g_param!=-999:
                        orbit_array = np.array([[0,semimaj_axis,eccentricity,np.deg2rad(inclination),np.deg2rad(long_ascnode),np.deg2rad(arg_perihelion),np.deg2rad(mean_anomaly),3,epoch_mjd,1,h_mag,g_param]],dtype=np.double,order='F')
                        mcd.output_log_entry(path_logfile,'Orbital elements successfully retrieved (ssobject_id={:d})'.format(ssobject_id))
                    else:
                        #orbit_array = np.array([[0,0,0,0,0,0,0,0,0,0,0,0]],dtype=np.double,order='F')
                        mcd.output_error_log_entry(path_logfile,path_errorfile,'Valid orbital elements not found for ssobject {:d}'.format(ssobject_id))
                        keep_going = False
                elif object_type == 'comet':
                    m1_mag          = 0               # temporary place-holder
                    k1_param        = 0               # temporary place-holder
                    perihelion_dist = float(row[5])
                    t_perihelion    = float(row[11])  ##### need to eventually convert t_perihelion to MJD #####
                    if perihelion_dist!=-999 and eccentricity!=-999 and inclination!=-999 and long_ascnode!=-999 and arg_perihelion!=-999 and t_perihelion!=-999 and epoch_mjd!=-999 and h_mag!=-999 and g_param!=-999:
                        orbit_array = np.array([[0,perihelion_dist,eccentricity,np.deg2rad(inclination),np.deg2rad(long_ascnode),np.deg2rad(arg_perihelion),t_perihelion,2,epoch_mjd,1,m1_mag,k1_param]],dtype=np.double,order='F')
                    else:
                        mcd.output_error_log_entry(path_logfile,path_errorfile,'Valid orbital elements not found for ssobject {:d}'.format(ssobject_id))
                        keep_going = False
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Valid object type ({:s}) not found for ssobject {:d}'.format(object_type,ssobject_id))
                    keep_going = False
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Object (ssobject_id={:d}) not found'.format(ssobject_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for ssobject {:d}: retrieve_orbital_elements()'.format(ssobject_id))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_orbital_elements()')
    return orbit_array,keep_going


def link_detection_temp_exposure(detection_id,search_result_date,search_result_telinst,search_result_datalink,exposures_not_found_file,sqlite_file,keep_going,path_logfile,path_errorfile):
    # Finds temporary exposure_id corresponding to specified detection for object position calculation
    #  (i.e., only retrieves first extension of full exposure to get observation time, without considering
    #   which extension contains the object or whether the object is even in any of the frames at all)
    temp_exposure_id,proc_data_path,temp_data_file = -1,'',''
    if keep_going:
        try:
            if get_filename_from_datalink(search_result_datalink)[6] != 'u':
                conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
                cursor1 = conn1.cursor()
                mcd.output_log_entry(path_logfile,'Searching for exposure corresponding to detection {:d}...'.format(detection_id))
                if search_result_telinst == 'SDSS':
                    raw_data_file = get_filename_from_datalink(search_result_datalink)[:-4]+'.fz'
                elif search_result_telinst == 'CFHT/MegaCam':
                    raw_data_file = get_filename_from_datalink(search_result_datalink)[:-3]+'.fz'
                query = "SELECT exposure_id,proc_data_path,proc_data_file FROM exposures WHERE raw_data_file = '{:s}'".format(raw_data_file)
                mcd.output_log_entry(path_logfile,query)
                cursor1.execute(query)
                row_exposure = cursor1.fetchone()
                if row_exposure != None:  # If corresponding exposure is found, identify specific extension with object
                    temp_exposure_id    = int(row_exposure[0])
                    proc_data_path      = row_exposure[1]
                    temp_data_file      = row_exposure[2]
                    mcd.output_log_entry(path_logfile,'Temporary exposure (exposure_id={:d}) found for detection_id {:d}'.format(temp_exposure_id,detection_id))
                else:
                    mcd.output_log_entry(path_logfile,'No exposures found for detection_id {:d}'.format(detection_id))
                    if os.path.exists(exposures_not_found_file):
                        with open(exposures_not_found_file,'a') as of:
                            of.write('{:s}\n'.format(get_filename_from_datalink(search_result_datalink)[:-4]+'.fz'))
                    else:
                        with open(exposures_not_found_file,'w') as of:
                            of.write('{:s}\n'.format(get_filename_from_datalink(search_result_datalink)[:-4]+'.fz'))
                    keep_going = False
                cursor1.close() # Close cursor
                conn1.close()   # Close connection to database file
            else:
                mcd.output_log_entry(path_logfile,'Skipping u-band file {:s}'.format(get_filename_from_datalink(search_result_datalink)[:-4]+'.fz'))
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for detection {:d}: link_detection_temp_exposure()'.format(detection_id))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: link_detection_temp_exposure()')
    return temp_exposure_id,proc_data_path,temp_data_file,keep_going
        
    
def compute_mjd_mid(exposure_id,sqlite_file,keep_going,path_logfile,path_errorfile):
    mjd_mid = -999
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            mcd.output_log_entry(path_logfile,'Retrieving MJD of midpoint of exposure {:d} from database...'.format(exposure_id))
            query = "SELECT exposure_start_jd,exposure_time FROM exposures WHERE exposure_id={:d}".format(exposure_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                exposure_start_jd = row[0]
                exposure_time     = row[1]
                mjd_mid = exposure_start_jd - 2400000.5 + (exposure_time/2/3600/24)
                keep_going = True
            else:
                mjd_mid = -999
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure information for Exposure {:d} not found'.format(exposure_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Exposure {:d}: compute_mjd_mid()'.format(exposure_id))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_mjd_mid()')
    return mjd_mid,keep_going


def compute_jd_mid(exposure_id,sqlite_file,keep_going,path_logfile,path_errorfile):
    jd_mid = -999
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            mcd.output_log_entry(path_logfile,'Retrieving JD of midpoint of exposure {:d} from database...'.format(exposure_id))
            query = "SELECT exposure_start_jd,exposure_time FROM exposures WHERE exposure_id={:d}".format(exposure_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                exposure_start_jd = row[0]
                exposure_time     = row[1]
                jd_mid = exposure_start_jd + (exposure_time/2/3600/24)
                keep_going = True
            else:
                jd_mid = -999
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure information for Exposure {:d} not found'.format(exposure_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for exposure {:d}: compute_jd_mid()'.format(exposure_id))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_jd_mid()')
    return jd_mid,keep_going


def retrieve_predicted_coords_horizons(desig_number,jd_mid,keep_going,path_logfile,path_errorfile):
    ra_deg,dec_deg = 0,0
    if keep_going:
        try:
            obs_code  = '500'   # observatory code
            mcd.output_log_entry(path_logfile,'Retrieving predicted RA/Dec for asteroid {:s} at {:.8f}'.format(desig_number,jd_mid))
            obj = Horizons(id=desig_number,location=obs_code,epochs=jd_mid)
            max_retries = 50
            for _ in range(max_retries):
                try:
                    ephems = obj.ephemerides()
                    print(obj.uri)
                    ra_deg      = ephems[0]['RA']
                    dec_deg     = ephems[0]['DEC']
                except:
                    mcd.output_log_entry(path_logfile,'JPL query failed. Retrying...')
                else:
                    break
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'JPL query failed for {:s},{:.8f}. Maximum retries reached'.format(desig_number,jd_mid))
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s},{:.8f}: retrieve_predicted_coords_horizons()'.format(desig_number,jd_mid))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_predicted_coords_horizons()')
    return ra_deg,dec_deg,keep_going

def identify_extension(proc_data_path,temp_data_file,right_ascension,declination,sqlite_file,keep_going,path_logfile,path_errorfile):
    extension_number_with_detection = -1
    if keep_going:
        try:
            fits_filestem = temp_data_file[:-16]  # e.g., frame-g-001045-5-0134_000_wcs.fits.fz --> frame-g-001045-5-0134
            os.chdir(proc_data_path)
            mcd.output_log_entry(path_logfile,'Processing {:s}{:s}'.format(proc_data_path,fits_filestem))
            for image_extension_to_check in glob.glob(fits_filestem+'_???_wcs.fits.fz'):
                extension_number = int(image_extension_to_check[:-12][-3:])  # get extension number from file name
                xpix,ypix = mcd.sky2pix_fzfile(image_extension_to_check,right_ascension,declination)
                if xpix != -1 and ypix != -1:
                    conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
                    cursor1 = conn1.cursor()
                    query = "SELECT exposure_npix_x,exposure_npix_y FROM exposures WHERE proc_data_path='{:s}' AND proc_data_file='{:s}'".format(proc_data_path,image_extension_to_check)
                    mcd.output_log_entry(path_logfile,query)
                    cursor1.execute(query)
                    row_exposure = cursor1.fetchone()
                    if row_exposure != None:  # If corresponding exposure is found, get image dimensions
                        npix_x = int(row_exposure[0])
                        npix_y = int(row_exposure[1])
                        mcd.output_log_entry(path_logfile,'Object at pixel coordinates {:f},{:f} in {:s}'.format(xpix,ypix,image_extension_to_check))
                        if xpix > 0 and xpix < npix_x and ypix > 0 and ypix < npix_y:
                            mcd.output_log_entry(path_logfile,'Detection found in extension {:d}'.format(extension_number))
                            extension_number_with_detection = extension_number
                        else:
                            mcd.output_log_entry(path_logfile,'Detection not in extension {:d} ({:d}x{:d})'.format(extension_number,npix_x,npix_y))
                    else:
                        mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure entry for {:s}{:s} not found'.format(image_extension_to_check))
                    cursor1.close() # Close cursor
                    conn1.close()   # Close connection to database file
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Unable to compute pixel coordinates for RA/Dec coordinates in {:s}'.format(image_extension_to_check))
            if extension_number_with_detection == -1:
                mcd.output_log_entry(path_logfile,'Object not found in any image extensions for {:s}.fits'.format(fits_filestem))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}.fits: identify_extension()'.format(fits_filestem))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: identify_extension()')
    return extension_number_with_detection,keep_going


def retrieve_exposure_id(extension_number_with_detection,proc_data_path,temp_data_file,sqlite_file,keep_going,path_logfile,path_errorfile):
    exposure_id = -1
    if keep_going:
        try:
            if extension_number_with_detection != -1:
                conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
                cursor1 = conn1.cursor()
                fits_filestem  = temp_data_file[:-16]  # e.g., frame-g-001045-5-0134_000_wcs.fits.fz --> frame-g-001045-5-0134
                proc_data_file = fits_filestem + '_{:03d}_wcs.fits.fz'.format(extension_number_with_detection)
                query = "SELECT exposure_id FROM exposures WHERE proc_data_path='{:s}' AND proc_data_file='{:s}'".format(proc_data_path,proc_data_file)
                mcd.output_log_entry(path_logfile,query)
                cursor1.execute(query)
                row_exposure = cursor1.fetchone()
                if row_exposure != None:  # If corresponding exposure is found, get exposure_id
                    exposure_id = row_exposure[0]
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure entry for {:s}{:s} not found'.format(proc_data_path,proc_data_file))
                    keep_going = False
                cursor1.close() # Close cursor
                conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}{:s}: retrieve_exposure_id()'.format(proc_data_path,proc_data_file))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_exposure_id()')
    return exposure_id,keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=6:
        print('Usage:\n python3 pyt_DT02_link_detection_data.py [base_path] [sqlite_file] [instrument] [thread_idx] [start_detection_id]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    instrument  = sys.argv[3]
    thread_idx          = int(sys.argv[4])
    start_detection_id  = int(sys.argv[5])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DT02_link_detection_data_{:s}_{:02d} execution started'.format(instrument,thread_idx),'DT02_link_detection_data_{:s}_{:02d} execution started.'.format(instrument,thread_idx))
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DT02_link_detection_data_{:s}_{:02d}'.format(instrument,thread_idx))

        # Create and open detection-exposure-object links and detection detail output files
        dir_detection_data  = base_path + 'detection_data/'
        if not os.path.isdir(dir_detection_data):
            mcd.create_directory(dir_detection_data,path_logfile,path_errorfile)
        output_file_links        = dir_detection_data + 'detection_exposure_object_links_{:s}_{:02d}_{:s}_toingest.txt'.format(instrument,thread_idx,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        exposures_not_found_file = dir_detection_data + 'exposures_not_found_{:s}_{:02d}_{:s}.txt'.format(instrument,thread_idx,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        with open(output_file_links,'w') as of:
            of.write('detection_id  exposure_id  ssobject_id\n')

        # Link detections with exposures and ssobjects
        link_detections_exposures_ssobjects(sqlite_file,output_file_links,exposures_not_found_file,instrument,thread_idx,start_detection_id,path_logfile,path_errorfile)

        # Compress output file
        mcd.compress_file_gzip(output_file_links)
        
    mcd.send_status_email('DT02_link_detection_data_{:s}_{:02d} execution complete'.format(instrument,thread_idx),'DT02_link_detection_data_{:s}_{:02d} execution complete.'.format(instrument,thread_idx))
    
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

