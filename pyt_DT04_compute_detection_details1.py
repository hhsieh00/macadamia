import re
import sys
import datetime
import math
import numpy as np
import glob, os, bz2, subprocess
import os.path
import sqlite3
from sqlite3 import Error
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.jplhorizons import Horizons
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DETECTION DETAIL COMPUTATION #####

def compute_detection_details(sqlite_file,output_file,thread_idx,start_detection_id,path_logfile,path_errorfile):
    try:
        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
        cursor = conn.cursor()
        
        mcd.output_log_entry(path_logfile,'Retrieving linked detection list...')
        query = "SELECT det.detection_id,dd.exposure_id,dd.ssobject_id,ss.desig_number FROM detections AS det JOIN detection_data AS dd JOIN ssobjects AS ss WHERE det.search_result_status='Exposure, object, and detection linked.' AND det.detection_id=dd.detection_id AND ss.ssobject_id=dd.ssobject_id ORDER BY det.detection_id ASC"
        mcd.output_log_entry(path_logfile,query)
        cursor.execute(query)
        rows = cursor.fetchall()
        mcd.output_log_entry(path_logfile,'{:d} linked detections retrieved'.format(len(rows)))
        for row in rows:
            keep_going   = True
            jd_mid_found = False
            detection_id = row[0]
            exposure_id  = row[1]
            ssobject_id  = row[2]
            desig_number = row[3]
            if '{:010d}'.format(detection_id)[-1:] == '{:01d}'.format(thread_idx) and detection_id > start_detection_id:
                mcd.output_log_entry(path_logfile,'Computing detection details for detection {:d}...'.format(detection_id))
            
                ### Compute mid-exposure MJD and retrieve orbital elements for specified object
                jd_mid,jd_mid_found,exposure_time,keep_going = compute_jd_mid(exposure_id,jd_mid_found,sqlite_file,keep_going,path_logfile,path_errorfile)

                ### Use Horizons to compute detection details
                if jd_mid_found:
                    keep_going = retrieve_write_detection_details_horizons(detection_id,desig_number,jd_mid,exposure_time,output_file,keep_going,path_logfile,path_errorfile)
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Mid-exposure JD not found for detection {:d}'.format(detection_id))
                    
        cursor.close() # Close cursor
        conn.close()   # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: compute_detection_details()')
        mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
        #print(e)
    return None


def compute_jd_mid(exposure_id,jd_mid_found,sqlite_file,keep_going,path_logfile,path_errorfile):
    # Compute mid-exposure MJD
    jd_mid = -999
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            mcd.output_log_entry(path_logfile,'Computing JD of midpoint of exposure {:d}...'.format(exposure_id))
            query = "SELECT exposure_start_jd,exposure_time FROM exposures WHERE exposure_id={:d}".format(exposure_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                exposure_start_jd = row[0]
                exposure_time     = row[1]
                jd_mid = exposure_start_jd + (exposure_time/2/3600/24)
                mcd.output_log_entry(path_logfile,'Mid-exposure JD for exposure {:d}: {:.8f}'.format(exposure_id,jd_mid))
                jd_mid_found = True
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure information for exposure_id {:d} not found'.format(exposure_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for exposure_id {:d}: compute_jd_mid()'.format(exposure_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            #print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_jd_mid()')
    return jd_mid,jd_mid_found,exposure_time,keep_going

def retrieve_write_detection_details_horizons(detection_id,desig_number,jd_mid,exposure_time,output_file,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Retrieving detection details for asteroid {:d}...'.format(detection_id))
            obs_code  = '500'   # observatory code
            obj = Horizons(id=desig_number,location=obs_code,epochs=jd_mid)
            mcd.output_log_entry_nonstring(path_logfile,obj)
            
            max_retries = 50
            for _ in range(max_retries):
                try:
                    ephems = obj.ephemerides()
                    #print(obj.uri)
                    mcd.output_log_entry_nonstring(path_logfile,obj.uri)
                    #print(ephems.columns)
                    ra_deg                  = ephems[0]['RA']
                    dec_deg                 = ephems[0]['DEC']
                    ra_err_deg              = ephems[0]['RA_3sigma']
                    dec_err_deg             = ephems[0]['DEC_3sigma']
                    ra_rate                 = ephems[0]['RA_rate']    # arcsec/hr
                    dec_rate                = ephems[0]['DEC_rate']   # arcsec/hr
                    tot_rate                = ((ra_rate**2)+(dec_rate**2))**0.5
                    trail_pa_expected       = np.arctan(ra_rate/dec_rate)/math.pi*180
                    trail_len_expected      = tot_rate / 3600 * exposure_time
                    phase_angle             = ephems[0]['alpha']
                    solar_elongation        = ephems[0]['elong']
                    heliocentric_dist       = ephems[0]['r']
                    geocentric_dist         = ephems[0]['delta']
                    mag_predicted           = ephems[0]['V']
                    pa_antisolar            = ephems[0]['sunTargetPA']
                    pa_neg_heliocentric_vel = ephems[0]['velocityPA']
                    ecliptic_longitude      = ephems[0]['EclLon']
                    ecliptic_latitude       = ephems[0]['EclLat']
                    lunar_illumination      = ephems[0]['lunar_illum']
                    lunar_elongation        = ephems[0]['lunar_elong']
                    true_anomaly            = ephems[0]['true_anom']
                    galactic_longitude      = ephems[0]['GlxLon']
                    galactic_latitude       = ephems[0]['GlxLat']
                    orbit_plane_angle       = ephems[0]['OrbPlaneAng']
                    with open(output_file,'a') as of:
                        #           det_id       ra        dec            ra_err         dec_err  ra_rate   dec_rate      trail_pa        trail_len   phase_angle     solar_elong    helio_dist  geo_dist          mag_pred     pa_antisolar           pa_neg_helio_vel  orbpl_angle       ecl_lon            ecl_lat            lunar_illum  lunar_elong       true_anomaly         gal_lon          gal_lat
                        of.write('  {:>10d}   {:11.7f}    {:11.7f}       {:11.7f}        {:11.7f}  {:13.7f}  {:13.7f}        {:7.2f}         {:7.2f}     {:8.4f}          {:8.4f}      {:13.8f}    {:13.8f}          {:5.2f}       {:7.3f}                  {:7.3f}     {:9.4f}            {:8.4f}           {:8.4f}              {:6.2f}            {:6.2f}      {:8.4f}            {:8.4f}           {:8.4f} \n'.format(detection_id,ra_deg,dec_deg,ra_err_deg,dec_err_deg,ra_rate,dec_rate,trail_pa_expected,trail_len_expected,phase_angle,solar_elongation,heliocentric_dist,geocentric_dist,mag_predicted,pa_antisolar,pa_neg_heliocentric_vel,orbit_plane_angle,ecliptic_longitude,ecliptic_latitude,lunar_illumination,lunar_elongation,true_anomaly,galactic_longitude,galactic_latitude))
                    mcd.output_log_entry(path_logfile,'Detection details for asteroid {:s} written to {:s}'.format(desig_number,output_file))
                    keep_going = True
                except:
                    mcd.output_log_entry(path_logfile,'JPL query failed for detection {:d}. Retrying...'.format(detection_id))
                else:
                    break
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'JPL query failed for detection {:d}'.format(detection_id))
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for detection {:d}: retrieve_write_detection_details_horizons()'.format(detection_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            #print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_write_detection_details_horizons()')
    return keep_going
    

def main():
    # Define filenames and paths
    if len(sys.argv)!=5:
        print('Usage:\n python3 pyt_ingest_detection_exposure_object_links.py [base_path] [sqlite_file] [thread_idx] [start_detection_id]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]
    thread_idx          = int(sys.argv[3])
    start_detection_id  = int(sys.argv[4])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DT04_compute_detection_details_{:01d} execution started'.format(thread_idx),'DT04_compute_detection_details_{:01d} execution started.'.format(thread_idx))
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'DT04_compute_detection_details_{:01d}'.format(thread_idx))

        dir_detection_data = base_path + 'detection_data/'
        if not os.path.isdir(dir_detection_data):
            mcd.create_directory(dir_detection_data,path_logfile,path_errorfile)

        output_file_details = dir_detection_data + 'detection_details_{:01d}_{:s}_toingest.txt'.format(thread_idx,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        with open(output_file_details,'w') as of:
            of.write('detection_id  ra_predicted  dec_predicted  ra_predicted_err  dec_predicted_err        ra_rate       dec_rate  trail_pa_pred  trail_len_pred  phase_angle  solar_elongation  heliocentric_dist  geocentric_dist  mag_predicted  pa_antisolar  pa_neg_heliocentric_vel  orb_pl_angle  ecliptic_longitude  ecliptic_latitude  lunar_illumination  lunar_elongation  true_anomaly  galactic_longitude  galactic_latitude\n')

        # Compute detection details
        compute_detection_details(sqlite_file,output_file_details,thread_idx,start_detection_id,path_logfile,path_errorfile)
        mcd.compress_file_gzip(output_file_details)

    mcd.send_status_email('DT04_compute_detection_details_{:01d} execution complete'.format(thread_idx),'DT04_compute_detection_details_{:01d} execution complete.'.format(thread_idx))
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

