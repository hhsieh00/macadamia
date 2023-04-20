import re
import sys
import time
import datetime
import math
import numpy as np
import glob, os, bz2, subprocess
import os.path
from numpy import linspace
from decimal import *
import sqlite3
from sqlite3 import Error
import macadamia_functions as mcd


def main():
    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_create_database.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = base_path + sys.argv[2]

    keep_going = True
    
    # Validate input parameters
    if base_path[-1:] != '/': base_path = base_path + '/'
    if not os.path.isdir(base_path):
        print('Directory {:s} not found.'.format(base_path))
        keep_going = False
    if os.path.isfile(sqlite_file):
        print('Database file {:s} already exists.'.format(sqlite_file))
        keep_going = False
    
    if keep_going:
        mcd.send_status_email('DB01_create_database execution started','{:s} - DB01_create_database execution started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        
        # Create and initialize log file
        path_logfile = mcd.initialize_log_file(base_path,'DB01_create_database')
            
        with open(path_logfile,'a') as log_file:

            # Connect to database file
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()

            print('{:s} - Creating asteroid photometry database...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - Creating asteroid photometry database...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

            # Create ssobjects table
            query = "CREATE TABLE IF NOT EXISTS ssobjects (ssobject_id INTEGER PRIMARY KEY, object_type TEXT NOT NULL, desig_number TEXT UNIQUE, desig_name TEXT, desig_provisional TEXT UNIQUE, h_mag REAL, g_param REAL, epoch_mjd REAL, semimaj_axis REAL, perihelion_dist REAL, eccentricity REAL, inclination REAL, arg_perihelion REAL, long_ascnode REAL, mean_anomaly REAL, t_perihelion REAL, date_added TEXT, date_updated TEXT)"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))

            # Create detections table
            query = "CREATE TABLE IF NOT EXISTS detections (detection_id INTEGER PRIMARY KEY, search_result_object TEXT, search_result_date TEXT, search_result_time TEXT, search_result_filter TEXT, search_result_exptime REAL, search_result_ra REAL, search_result_dec REAL, search_result_target TEXT, search_result_telinst TEXT, search_result_datalink TEXT, search_result_status TEXT, CONSTRAINT search_result UNIQUE (search_result_object,search_result_date,search_result_time,search_result_filter,search_result_telinst))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            # Create detection_data table
            query = "CREATE TABLE IF NOT EXISTS detection_data (detection_id INTEGER UNIQUE, exposure_id INTEGER, ssobject_id INTEGER, x_pixcoord_predicted REAL, y_pixcoord_predicted REAL, x_pixcoord REAL, y_pixcoord REAL, right_ascension REAL, declination REAL, ra_err REAL, dec_err REAL, dist_predicted_position REAL, flux REAL, dflux REAL, mag REAL, mag_err REAL, psf_fwhm REAL, signal_to_noise REAL, trail_length REAL, trail_fwhm REAL, trail_phi REAL, ra_source_1 REAL, dec_source_1 REAL, dist_source_pred_1 REAL, dist_source_actual_1 REAL, mag_source_1 REAL, magerr_source_1 REAL, ra_source_2 REAL, dec_source_2 REAL, dist_source_pred_2 REAL, dist_source_actual_2 REAL, mag_source_2 REAL, magerr_source_2 REAL, ra_source_3 REAL, dec_source_3 REAL, dist_source_pred_3 REAL, dist_source_actual_3 REAL, mag_source_3 REAL, magerr_source_3 REAL, source_density_1arcmin REAL, dist_edge_left REAL, dist_edge_right REAL, dist_edge_bottom REAL, dist_edge_top REAL, preview_image_path TEXT, preview_image_file TEXT, detection_status TEXT, FOREIGN KEY (detection_id) REFERENCES detections(detection_id), FOREIGN KEY (exposure_id) REFERENCES exposures(exposure_id), FOREIGN KEY (ssobject_id) REFERENCES ssobjects(ssobject_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            # Create detection_details table
            query = "CREATE TABLE IF NOT EXISTS detection_details (detection_id INTEGER UNIQUE, ra_predicted REAL, dec_predicted REAL, ra_predicted_err REAL, dec_predicted_err REAL, ra_rate REAL, dec_rate REAL, trail_pa_expected REAL, trail_length_expected REAL, ecliptic_longitude REAL, ecliptic_latitude REAL, heliocentric_dist REAL, geocentric_dist REAL, solar_elongation REAL, phase_angle REAL, lunar_elongation REAL, lunar_illumination REAL, pa_antisolar REAL, pa_neg_heliocentric_vel REAL, orb_plane_angle REAL, galactic_longitude REAL, galactic_latitude REAL, true_anomaly REAL, vmag_predicted REAL, vmag_predicted_err REAL, detection_details_status TEXT, FOREIGN KEY (detection_id) REFERENCES detections(detection_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            # Create detection_multiap_photometry table
            query = "CREATE TABLE IF NOT EXISTS detection_multiap_photometry (detection_id INTEGER UNIQUE, flux_t_05rKron REAL, dflux_t_05rKron REAL, mag_t_05rKron REAL, magerr_t_05rKron REAL, flux_t_10rKron REAL, dflux_t_10rKron REAL, mag_t_10rKron REAL, magerr_t_10rKron REAL, flux_t_20rKron REAL, dflux_t_20rKron REAL, mag_t_20rKron REAL, magerr_t_20rKron REAL, flux_t_30rKron REAL, dflux_t_30rKron REAL, mag_t_30rKron REAL, magerr_t_30rKron REAL, flux_t_40rKron REAL, dflux_t_40rKron REAL, mag_t_40rKron REAL, magerr_t_40rKron REAL, flux_t_50rKron REAL, dflux_t_50rKron REAL, mag_t_50rKron REAL, magerr_t_50rKron REAL, FOREIGN KEY (detection_id) REFERENCES detections(detection_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))

            # Create exposures table
            query = "CREATE TABLE IF NOT EXISTS exposures (exposure_id INTEGER PRIMARY KEY, mosaic_element_id INTEGER, date_tai TEXT NOT NULL, exposure_start_tai TEXT NOT NULL, exposure_start_jd REAL, exposure_time REAL, exposure_npix_x INTEGER, exposure_npix_y INTEGER, filter_name TEXT, airmass REAL, pointing_center_ra REAL, pointing_center_dec REAL, tracking_rate_ra REAL, tracking_rate_dec REAL, binning_x INTEGER, binning_y INTEGER, wcs_nstars INTEGER, wcs_crpix1 REAL, wcs_crpix2 REAL, wcs_crval1 REAL, wcs_crval2 REAL, wcs_cd1_1 REAL, wcs_cd1_2 REAL, wcs_cd2_1 REAL, wcs_cd2_2 REAL, wcs_pixscale_x REAL, wcs_pixscale_y REAL, sky_mean REAL, sky_stddev REAL, zero_point REAL, zpoint_err REAL, zpoint_nstars INTEGER, limiting_mag_ps_3sig REAL, limiting_mag_sb_3sig REAL, psf_width_mean REAL, source_density REAL, base_filename TEXT, raw_data_link TEXT, raw_data_path TEXT, raw_data_file TEXT, proc_data_path TEXT, proc_data_file, source_table_path, source_table_file, exposure_status TEXT, FOREIGN KEY (mosaic_element_id) REFERENCES mosaic_elements(mosaic_element_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            # Create mosaic_elements table
            query = "CREATE TABLE IF NOT EXISTS mosaic_elements (mosaic_element_id INTEGER PRIMARY KEY, instrument_id INTEGER, mosaic_element_num INTEGER, gain REAL, read_noise REAL, npix_x INTEGER, npix_y INTEGER, pixel_scale_x REAL, pixel_scale_y REAL, FOREIGN KEY (instrument_id) REFERENCES instruments(instrument_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            # Create instruments table
            query = "CREATE TABLE IF NOT EXISTS instruments (instrument_id INTEGER PRIMARY KEY, telescope_id INTEGER, instrument_name TEXT, FOREIGN KEY (telescope_id) REFERENCES telescopes(telescope_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))

            # Create telescopes table
            query = "CREATE TABLE IF NOT EXISTS telescopes (telescope_id INTEGER PRIMARY KEY, telescope_name TEXT, observatory_code TEXT, aperture_size REAL, latitude REAL, longitude REAL, elevation REAL)"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))

            # Create search_history table
            query = "CREATE TABLE IF NOT EXISTS search_history (ssobject_id INTEGER, instrument_id INTEGER, last_searched TEXT, FOREIGN KEY (ssobject_id) REFERENCES ssobjects(ssobject_id), FOREIGN KEY (instrument_id) REFERENCES instruments(instrument_id))"
            cursor.execute(query)
            log_file.write("{:s} - {:s}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),query))
        
            print('{:s} - >>> Asteroid photometry database created.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Asteroid photometry database created.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
            # Committing changes and closing the connection to the database file
            conn.commit()
            conn.close()
            print('{:s} - >>> Database closed.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Database closed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

        mcd.send_status_email('DB01_create_database execution complete','{:s} - DB01_create_database execution complete.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None

        
if __name__ == '__main__':
    main()

