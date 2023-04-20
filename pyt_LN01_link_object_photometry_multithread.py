import re
import sys
import time
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
import glob, os, bz2, subprocess
import os.path
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from decimal import *
from scipy.stats import norm
import sqlite3
from sqlite3 import Error
import jdcal
from jdcal import gcal2jd,jd2gcal
import statistics
from uncertainties import unumpy,ufloat
from uncertainties.umath import *
from uncertainties.unumpy import *
import montage_wrapper as montage
from astropy.nddata.utils import Cutout2D
import astropy.visualization as vis
import astropy.units as un
from astroquery.jplhorizons import Horizons
import macadamia_functions as mcd

### Edit History
# 2021-07-01: implemented new logging functionality
# 2021-09-23: changed initial SQL query to include already processed images to update database/previews;
#              changed error handling to print "e" to logfile


##### FUNCTION DEFINITIONS -- DETECTION DETAIL COMPUTATION #####

def link_object_photometry_multithread(thread_idx,start_detection_id,sqlite_file,output_data_file,preview_path_stem,preview_size_arcsec,path_logfile,path_errorfile):
    try:
        conn = mcd.create_connection(sqlite_file)  # Open connection to database file
        cursor = conn.cursor()
        mcd.output_log_entry(path_logfile,'Retrieving linked detection list...')
        query = "SELECT detection_id,exposure_id FROM detection_data WHERE detection_status='Exposure, object, and detection linked.' ORDER BY detection_id"
        mcd.output_log_entry(path_logfile,query)
        cursor.execute(query)
        rows = cursor.fetchall()
        mcd.output_log_entry(path_logfile,'{:d} linked detections retrieved'.format(len(rows)))
        for row in rows:
            keep_going = True
            detection_id = row[0]
            exposure_id  = row[1]
            try:
                if '{:010d}'.format(detection_id)[-2:] == '{:02d}'.format(thread_idx) and detection_id > start_detection_id:
                    mcd.output_log_entry(path_logfile,'Linking object photometry for detection {:d}...'.format(detection_id))
                    proc_data_fullpath_fz,source_table_fullpath,preview_size_pix,pixscale_x,npix_x,npix_y,keep_going = retrieve_exposure_data(exposure_id,preview_size_arcsec,sqlite_file,keep_going,path_logfile,path_errorfile)
                    ra_predicted,dec_predicted,keep_going = retrieve_predicted_radec(detection_id,sqlite_file,keep_going,path_logfile,path_errorfile)
                    field_source_table,keep_going = create_field_source_table(source_table_fullpath,keep_going,path_logfile,path_errorfile)
                    field_source_table,target_found,preview_path,preview_filename,keep_going = find_target_source(detection_id,ra_predicted,dec_predicted,pixscale_x,field_source_table,proc_data_fullpath_fz,preview_size_pix,preview_path_stem,keep_going,path_logfile,path_errorfile)
                    keep_going = compute_write_detection_data(detection_id,proc_data_fullpath_fz,ra_predicted,dec_predicted,pixscale_x,npix_x,npix_y,field_source_table,target_found,output_data_file,preview_path,preview_filename,keep_going,path_logfile,path_errorfile)
                    mcd.output_log_entry(path_logfile,'Linking object photometry done.')
            except Exception as e:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Detection {:d}, Exposure {:d}: link_object_photometry_multithread()'.format(detection_id,exposure_id))
                mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
        cursor.close() # Close cursor
        conn.close()   # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: link_object_photometry_multithread()')
        mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
    return None


def compute_write_detection_data(detection_id,proc_data_fullpath_fz,ra_predicted,dec_predicted,pixscale_x,npix_x,npix_y,field_source_table,target_found,output_data_file,preview_path,preview_filename,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            dist_source_actual_1 = -999
            dist_source_actual_2 = -999
            dist_source_actual_3 = -999
            mcd.output_log_entry(path_logfile,'Writing detection data to detection_data table for Detection {:d}...'.format(detection_id))
            num_sources = len(field_source_table)
            x_pixcoord_predicted,y_pixcoord_predicted = mcd.sky2pix_fzfile(proc_data_fullpath_fz,ra_predicted,dec_predicted)
            if target_found:
                x_pixcoord              = field_source_table[0]['x']
                y_pixcoord              = field_source_table[0]['y']
                right_ascension         = field_source_table[0]['RA']
                declination             = field_source_table[0]['Dec']
                dist_predicted_position = field_source_table[0]['target_dist']
                flux_s                  = field_source_table[0]['flux_s']
                dflux_s                 = field_source_table[0]['dflux_s']
                flux_t                  = field_source_table[0]['flux_t']
                dflux_t                 = field_source_table[0]['dflux_t']
                mag_s                   = field_source_table[0]['mag_s']
                magerr_s                = field_source_table[0]['magerr_s']
                mag_t                   = field_source_table[0]['mag_t']
                magerr_t                = field_source_table[0]['magerr_t']
                SNR_t                   = field_source_table[0]['SNR_t']
                trail_t                 = field_source_table[0]['trail_t']*pixscale_x
                fwhm_t                  = field_source_table[0]['fwhm_t']*pixscale_x
                fwhm_s                  = field_source_table[0]['fwhm_s']*pixscale_x
                phi_t                   = field_source_table[0]['phi_t']
                flux_t_05rK             = field_source_table[0]['flux_t_05rK']
                dflux_t_05rK            = field_source_table[0]['dflux_t_05rK']
                mag_t_05rK              = field_source_table[0]['mag_t_05rK']
                magerr_t_05rK           = field_source_table[0]['magerr_t_05rK']
                flux_t_10rK             = field_source_table[0]['flux_t_10rK']
                dflux_t_10rK            = field_source_table[0]['dflux_t_10rK']
                mag_t_10rK              = field_source_table[0]['mag_t_10rK']
                magerr_t_10rK           = field_source_table[0]['magerr_t_10rK']
                flux_t_20rK             = field_source_table[0]['flux_t_20rK']
                dflux_t_20rK            = field_source_table[0]['dflux_t_20rK']
                mag_t_20rK              = field_source_table[0]['mag_t_20rK']
                magerr_t_20rK           = field_source_table[0]['magerr_t_20rK']
                flux_t_30rK             = field_source_table[0]['flux_t_30rK']
                dflux_t_30rK            = field_source_table[0]['dflux_t_30rK']
                mag_t_30rK              = field_source_table[0]['mag_t_30rK']
                magerr_t_30rK           = field_source_table[0]['magerr_t_30rK']
                flux_t_40rK             = field_source_table[0]['flux_t_40rK']
                dflux_t_40rK            = field_source_table[0]['dflux_t_40rK']
                mag_t_40rK              = field_source_table[0]['mag_t_40rK']
                magerr_t_40rK           = field_source_table[0]['magerr_t_40rK']
                flux_t_50rK             = field_source_table[0]['flux_t_50rK']
                dflux_t_50rK            = field_source_table[0]['dflux_t_50rK']
                mag_t_50rK              = field_source_table[0]['mag_t_50rK']
                magerr_t_50rK           = field_source_table[0]['magerr_t_50rK']
                if num_sources > 1:
                    dist_source_actual_1 = compute_angular_dist(right_ascension,declination,field_source_table[1]['RA'],field_source_table[1]['Dec'])*3600
                else:
                    dist_source_actual_1 = -999
                if num_sources > 2:
                    dist_source_actual_2 = compute_angular_dist(right_ascension,declination,field_source_table[2]['RA'],field_source_table[2]['Dec'])*3600
                else:
                    dist_source_actual_2 = -999
                if num_sources > 3:
                    dist_source_actual_3 = compute_angular_dist(right_ascension,declination,field_source_table[3]['RA'],field_source_table[3]['Dec'])*3600
                else:
                    dist_source_actual_3 = -999
                idx = 1
            else:
                x_pixcoord,y_pixcoord       = -999,-999
                right_ascension,declination = -999,-999
                dist_predicted_position     = -999
                flux_s,dflux_s              = -999,-999
                flux_t,dflux_t              = -999,-999
                mag_s,magerr_s              = 99.99,99.99
                mag_t,magerr_t              = 99.99,99.99
                SNR_t                       = -999
                trail_t                     = -999
                fwhm_t,fwhm_s               = -999,-999
                phi_t                       = -999
                flux_t_05rK,dflux_t_05rK    = -999,-999
                flux_t_10rK,dflux_t_10rK    = -999,-999
                flux_t_20rK,dflux_t_20rK    = -999,-999
                flux_t_30rK,dflux_t_30rK    = -999,-999
                flux_t_40rK,dflux_t_40rK    = -999,-999
                flux_t_50rK,dflux_t_50rK    = -999,-999
                mag_t_05rK,magerr_t_05rK    = 99.99,99.99
                mag_t_10rK,magerr_t_10rK    = 99.99,99.99
                mag_t_20rK,magerr_t_20rK    = 99.99,99.99
                mag_t_30rK,magerr_t_30rK    = 99.99,99.99
                mag_t_40rK,magerr_t_40rK    = 99.99,99.99
                mag_t_50rK,magerr_t_50rK    = 99.99,99.99
                idx                     = 0
            if num_sources > idx:
                ra_source_1        = field_source_table[idx]['RA']
                dec_source_1       = field_source_table[idx]['Dec']
                dist_source_pred_1 = field_source_table[idx]['target_dist']
                mag_source_1       = field_source_table[idx]['mag_s']
                magerr_source_1    = field_source_table[idx]['magerr_s']
            else:
                ra_source_1        = -999
                dec_source_1       = -999
                dist_source_pred_1 = -999
                mag_source_1       = 99.99
                magerr_source_1    = 99.99
            if num_sources > (idx+1):
                ra_source_2        = field_source_table[idx+1]['RA']
                dec_source_2       = field_source_table[idx+1]['Dec']
                dist_source_pred_2 = field_source_table[idx+1]['target_dist']
                mag_source_2       = field_source_table[idx+1]['mag_s']
                magerr_source_2    = field_source_table[idx+1]['magerr_s']
            else:
                ra_source_2        = -999
                dec_source_2       = -999
                dist_source_pred_2 = -999
                mag_source_2       = 99.99
                magerr_source_2    = 99.99
            if num_sources > (idx+2):
                ra_source_3        = field_source_table[idx+2]['RA']
                dec_source_3       = field_source_table[idx+2]['Dec']
                dist_source_pred_3 = field_source_table[idx+2]['target_dist']
                mag_source_3       = field_source_table[idx+2]['mag_s']
                magerr_source_3    = field_source_table[idx+2]['magerr_s']
            else:
                ra_source_3        = -999
                dec_source_3       = -999
                dist_source_pred_3 = -999
                mag_source_3       = 99.99
                magerr_source_3    = 99.99
            num_sources_1arcmin = 0
            while field_source_table[idx]['target_dist'] < 60:
                num_sources_1arcmin += 1
                idx += 1
            source_density_1arcmin = num_sources_1arcmin / (1*1*math.pi)  # density of sources as measured within r=1arcmin circle
            dist_edge_left   = x_pixcoord_predicted * pixscale_x
            dist_edge_right  = npix_x - x_pixcoord_predicted * pixscale_x
            dist_edge_bottom = y_pixcoord_predicted * pixscale_x
            dist_edge_top    = npix_y - y_pixcoord_predicted * pixscale_x
            with open(output_data_file,'a') as of:
                #         detid       x_pred    y_pred   x_coord  y_coord    ra        dec               dist_p   flux_s    dflux_s   flux_t   dflux_t    mag_s  magerr_s   mag_t  magerr_t   fwhm_s    SNR       trail     fwhm_t    phi     f_05rk    df_05rk       m_05rk          me_05rK  f_10rk    df_10rk       m_10rk          me_10rK  f_20rk    df_20rk       m_20rk          me_20rK  f_30rk    df_30rk       m_30rk          me_30rK  f_40rk    df_40rk       m_40rk          me_40rK  f_50rk    df_50rk       m_50rk          me_50rK    ra1       dec1    dist1_p   dist1_a      mag1       magerr1    ra2       dec2    dist2_p   dist2_a      mag2       magerr2    ra3       dec3    dist3_p   dist3_a      mag3       magerr3         srcdensity   dist_l   dist_r    dist_b    dist_t    file  dir
                of.write('{:010d}    {:12.3f}  {:12.3f}  {:9.3f}  {:9.3f}    {:13.8f}  {:13.8f}         {:8.3f}  {:12.6f}  {:12.6f}  {:12.6f}  {:12.6f}  {:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}  {:9.3f}   {:8.3f}    {:8.3f}   {:8.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:11.3f}  {:12.3f}      {:6.3f}         {:6.3f}  {:13.8f}  {:13.8f}  {:15.3f}  {:17.3f}   {:8.3f}      {:8.3f}  {:13.8f}  {:13.8f}  {:15.3f}  {:17.3f}   {:8.3f}      {:8.3f}  {:13.8f}  {:13.8f}  {:15.3f}  {:17.3f}   {:8.3f}      {:8.3f}           {:10.3f}  {:14.1f}  {:15.1f}  {:16.1f}  {:13.1f}  {:s}  {:s}\n'.format(detection_id,x_pixcoord_predicted,y_pixcoord_predicted,x_pixcoord,y_pixcoord,right_ascension,declination,dist_predicted_position,flux_s,dflux_s,flux_t,dflux_t,mag_s,magerr_s,mag_t,magerr_t,fwhm_s,SNR_t,trail_t,fwhm_t,phi_t,flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK,flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK,flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK,flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK,flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK,flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK,ra_source_1,dec_source_1,dist_source_pred_1,dist_source_actual_1,mag_source_1,magerr_source_1,ra_source_2,dec_source_2,dist_source_pred_2,dist_source_actual_2,mag_source_2,magerr_source_2,ra_source_3,dec_source_3,dist_source_pred_3,dist_source_actual_3,mag_source_3,magerr_source_3,source_density_1arcmin,dist_edge_left,dist_edge_right,dist_edge_bottom,dist_edge_top,preview_filename+'.png',preview_path))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Detection {:d}: compute_write_detection_data()'.format(detection_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_write_detection_data()')
    return keep_going


def find_target_source(detection_id,ra_predicted,dec_predicted,pixscale_x,field_source_table,proc_data_fullpath_fz,preview_size_pix,preview_path_stem,keep_going,path_logfile,path_errorfile):
    preview_path,preview_filename = '',''
    target_found = False
    if keep_going:
        try:
            num_sources = len(field_source_table)
            # Identify nearest sources to target object's predicted position
            for idx in range(0,num_sources):  # compute distance in arcsec from each source to target object's predicted position
                field_source_table[idx]['target_dist'] = 3600 * compute_angular_dist(field_source_table[idx]['RA'],field_source_table[idx]['Dec'],ra_predicted,dec_predicted)
            field_source_table.sort('target_dist')  # sort calibration star table in increasing order of distance to target
            if field_source_table[0]['target_dist'] < 5:
                mcd.output_log_entry(path_logfile,'Target detected at {:.1f},{:.1f} ({:.1f} arcsec from predicted position)'.format(field_source_table[0]['x'],field_source_table[0]['y'],field_source_table[0]['target_dist']))
                mcd.output_log_entry(path_logfile,'Target:  x={:.1f}; y={:.1f}; target_dist={:.1f}, flux={:.1f}; SNR={:.1f}; fwhm_width={:.2f} arcsec ({:.2f} pix); mag={:.3f}+/-{:.3f}'.format(field_source_table[0]['x'],field_source_table[0]['y'],field_source_table[0]['target_dist'],field_source_table[0]['flux_t'],field_source_table[0]['SNR_t'],field_source_table[0]['fwhm_t']*pixscale_x,field_source_table[0]['fwhm_t'],field_source_table[0]['mag_t'],field_source_table[0]['magerr_t']))
                preview_path,preview_filename,keep_going = make_preview_image(proc_data_fullpath_fz,preview_path_stem,detection_id,field_source_table[0]['RA'],field_source_table[0]['Dec'],preview_size_pix,keep_going,path_logfile,path_errorfile)
                target_found = True
                keep_going = True
            else:
                mcd.output_log_entry(path_logfile,'Target object not detected for Detection {:d}'.format(detection_id))
                preview_path,preview_filename,keep_going = make_preview_image(proc_data_fullpath_fz,preview_path_stem,detection_id,ra_predicted,dec_predicted,preview_size_pix,keep_going,path_logfile,path_errorfile)
                target_found = False
                keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Detection {:d}: find_target_source()'.format(detection_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: find_target_source()')
    return field_source_table,target_found,preview_path,preview_filename,keep_going

    
def retrieve_exposure_data(exposure_id,preview_size_arcsec,sqlite_file,keep_going,path_logfile,path_errorfile):
    proc_data_fullpath_fz,source_table_fullpath = '',''
    preview_size_pix,wcs_pixscale_x,exposure_npix_x,exposure_npix_y = -1,0,0,0
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            query = "SELECT proc_data_path,proc_data_file,source_table_path,source_table_file,wcs_pixscale_x,exposure_npix_x,exposure_npix_y,exposure_status FROM exposures WHERE exposure_id={:d}".format(exposure_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                proc_data_path        = row[0]
                proc_data_file        = row[1]
                source_table_path     = row[2]
                source_table_file     = row[3]
                wcs_pixscale_x        = float(row[4])
                npix_x                = int(row[5])
                npix_y                = int(row[6])
                exposure_status       = row[7]
                mcd.output_log_entry(path_logfile,'Data for exposure {:d} found'.format(exposure_id))
                if not (exposure_status == 'Photometric calibration successful.' or exposure_status == 'Photometric calibration successful' or source_table_file != None):
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Photometric calibration of exposure {:d} not completed successfully or not yet ingested'.format(exposure_id))
                    keep_going = False
                if proc_data_path != None and proc_data_file != None:
                    proc_data_fullpath_fz = proc_data_path + proc_data_file
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Processed data path not found for exposure {:d}'.format(exposure_id))
                    keep_going = False
                if source_table_path != None and source_table_file != None:
                    source_table_fullpath = source_table_path + source_table_file
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Calibrated source table not found for exposure {:d}'.format(exposure_id))
                    keep_going = False
                if preview_size_arcsec != None and wcs_pixscale_x != None and npix_x != None and npix_y != None:
                    preview_size_pix = preview_size_arcsec / wcs_pixscale_x
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure parameters not found for exposure {:d}'.format(exposure_id))
                    keep_going = False
                if exposure_status == None and source_table_file == None:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure status not found for exposure {:d}'.format(exposure_id))
                    keep_going = False
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Exposure {:d} not found'.format(exposure_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for exposure {:d}: retrieve_exposure_data()'.format(exposure_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_exposure_data()')
    return proc_data_fullpath_fz,source_table_fullpath,preview_size_pix,wcs_pixscale_x,npix_x,npix_y,keep_going


def retrieve_predicted_radec(detection_id,sqlite_file,keep_going,path_logfile,path_errorfile):
    ra_predicted,dec_predicted = 0,0
    if keep_going:
        try:
            conn1 = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor1 = conn1.cursor()
            query = "SELECT ra_predicted,dec_predicted FROM detection_details WHERE detection_id={:d}".format(detection_id)
            mcd.output_log_entry(path_logfile,query)
            cursor1.execute(query)
            row = cursor1.fetchone()
            if row != None:
                ra_predicted  = row[0]
                dec_predicted = row[1]
                mcd.output_log_entry(path_logfile,'Predicted RA/Dec coordinates for detection {:d} found.'.format(detection_id))
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Detection {:d} not found'.format(detection_id))
                keep_going = False
            cursor1.close() # Close cursor
            conn1.close()   # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Detection {:d}: retrieve_predicted_radec()'.format(detection_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_predicted_radec()')
    return ra_predicted,dec_predicted,keep_going


def create_field_source_table(source_table_fullpath,keep_going,path_logfile,path_errorfile):
    # Create calibration star table
    field_source_table = Table()
    if keep_going:
        try:
            multiap_photom_data_fullpath = source_table_fullpath[:-23]+'_multiap_photometry.txt'
            mcd.output_log_entry(path_logfile,'Creating field source table from {:s} and {:s}...'.format(source_table_fullpath,multiap_photom_data_fullpath))
            source_table_fullpath_gz        = source_table_fullpath + '.gz'
            multiap_photom_data_fullpath_gz = multiap_photom_data_fullpath + '.gz'
            if (os.path.exists(source_table_fullpath_gz) or os.path.exists(source_table_fullpath)) and (os.path.exists(multiap_photom_data_fullpath_gz) or os.path.exists(multiap_photom_data_fullpath)):
                if os.path.exists(source_table_fullpath_gz) and not os.path.exists(source_table_fullpath):
                    mcd.decompress_file_gzip(source_table_fullpath_gz)
                if os.path.exists(multiap_photom_data_fullpath_gz) and not os.path.exists(multiap_photom_data_fullpath):
                    mcd.decompress_file_gzip(multiap_photom_data_fullpath_gz)
                x                  = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(0))
                y                  = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(1))
                flux_s             = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(2))
                dflux_s            = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(3))
                major_s            = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(4))
                minor_s            = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(5))
                flux_t             = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(11))
                dflux_t            = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(12))
                trail_t            = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(13))
                fwhm_t             = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(14))
                phi_t              = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(15))
                chiN_t             = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(16))
                SNR_t              = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(17))
                RA                 = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(18))
                Dec                = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(19))
                mag_s              = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(21))
                magerr_s           = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(22))
                mag_t              = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(23))
                magerr_t           = np.genfromtxt(source_table_fullpath,skip_header=1,usecols=(24))
                flux_t_05rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(4))
                dflux_t_05rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(5))
                mag_t_05rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(6))
                magerr_t_05rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(7))
                flux_t_10rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(10))
                dflux_t_10rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(11))
                mag_t_10rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(12))
                magerr_t_10rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(13))
                flux_t_20rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(16))
                dflux_t_20rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(17))
                mag_t_20rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(18))
                magerr_t_20rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(19))
                flux_t_30rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(22))
                dflux_t_30rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(23))
                mag_t_30rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(24))
                magerr_t_30rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(25))
                flux_t_40rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(28))
                dflux_t_40rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(29))
                mag_t_40rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(30))
                magerr_t_40rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(31))
                flux_t_50rK        = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(34))
                dflux_t_50rK       = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(35))
                mag_t_50rK         = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(36))
                magerr_t_50rK      = np.genfromtxt(multiap_photom_data_fullpath,skip_header=1,usecols=(37))
                numpts = len(major_s)
                fwhm_s = [0 for idx in range(numpts)]
                for idx in range(0,numpts):
                    fwhm_s[idx] = ((major_s[idx]**2 + minor_s[idx]**2)/2)**0.5
                field_source_table = Table([x,y,flux_s,dflux_s,flux_t,dflux_t,trail_t,fwhm_t,phi_t,chiN_t,SNR_t,fwhm_s,RA,Dec,mag_s,magerr_s,mag_t,magerr_t,flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK,flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK,flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK,flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK,flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK,flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK], \
                    names=('x','y','flux_s','dflux_s','flux_t','dflux_t','trail_t','fwhm_t','phi_t','chiN_t','SNR_t','fwhm_s','RA','Dec','mag_s','magerr_s','mag_t','magerr_t','flux_t_05rK','dflux_t_05rK','mag_t_05rK','magerr_t_05rK','flux_t_10rK','dflux_t_10rK','mag_t_10rK','magerr_t_10rK','flux_t_20rK','dflux_t_20rK','mag_t_20rK','magerr_t_20rK','flux_t_30rK','dflux_t_30rK','mag_t_30rK','magerr_t_30rK','flux_t_40rK','dflux_t_40rK','mag_t_40rK','magerr_t_40rK','flux_t_50rK','dflux_t_50rK','mag_t_50rK','magerr_t_50rK'), \
                    dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
                field_source_table['target_dist'] = -999.0
                mcd.output_log_entry(path_logfile,'Creating field source table done.')
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Field source table files ({:s} and {:s}) not found'.format(source_table_fullpath,multiap_photom_data_fullpath))
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: create_field_source_table()'.format(source_table_fullpath))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_field_source_table()')
    return field_source_table,keep_going
            

def compute_angular_dist(ra1,dec1,ra2,dec2):
    # Compute angular distance in degrees from two sets of RA,Dec coordinates, also in degrees
    source_distance_deg = 0
    try:
        c1 = SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
        c2 = SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
        source_distance_deg = c1.separation(c2).degree
    except Exception as e:
        print('Function failed: compute_angular_dist()')
        print(e)
    return source_distance_deg


def make_preview_image(imagename,preview_path_stem,detection_id,target_ra,target_dec,preview_size_pix,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Creating preview image for Detection {:d}'.format(detection_id))
            mcd.decompress_file_funpack(imagename)
            imagename = imagename[:-3]
            # generate preview image path
            preview_dir_lowerbound = (floor((detection_id-1) / 10000))*10000 + 1
            preview_dir_upperbound = preview_dir_lowerbound + 9999
            preview_path = preview_path_stem + '_{:010d}_{:010d}/'.format(preview_dir_lowerbound,preview_dir_upperbound)
            mcd.create_directory(preview_path,path_logfile,path_errorfile)
            preview_filename = 'preview_det{:010d}'.format(detection_id) # set output filename
            image,header,wcs = open_image(imagename)                                    # load image and align with North up
            target_xcoord,target_ycoord = wcs.wcs_world2pix(target_ra,target_dec,1)     # find pixel coords of target in new image
    
            # Check image orientation
            test_ra1,test_dec1 = wcs.wcs_pix2world(1,1,1)
            test_ra2,test_dec2 = wcs.wcs_pix2world(2,1,1)
            if test_ra1 > test_dec2:
                flip = False
            else:
                flip = True
    
            # Select a region centered on the target
            cutout = Cutout2D(image,(target_xcoord,target_ycoord),(preview_size_pix,preview_size_pix),mode='partial',fill_value=0)

            # Set image scale by measuring stats of full image (not just cutout)
            imagedata_f = image.flatten()
            pix_median = statistics.median(imagedata_f)
            mu,sigma = norm.fit(imagedata_f)

            if flip:
                plt.imshow(cutout.data[:,::-1],aspect='equal',origin='lower',cmap='gray',norm=mcolors.PowerNorm(gamma=0.9,vmin=np.amin(imagedata_f),vmax=(pix_median+0.5*sigma)))
            else:
                plt.imshow(cutout.data,aspect='equal',origin='lower',cmap='gray',norm=mcolors.PowerNorm(gamma=0.9,vmin=np.amin(imagedata_f),vmax=(pix_median+0.5*sigma)))
            plt.axis('off')
            plt.savefig(preview_path + preview_filename + '.png',format='png',transparent=True)
            plt.show()
            mcd.output_log_entry(path_logfile,'Saving preview image to {:s}'.format(preview_path + preview_filename+'.png'))
            plt.clf()
            plt.cla()
            plt.close()
            mcd.compress_file_fpack(imagename)
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for Detection {:d}: make_preview_image()'.format(detection_id))
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: make_preview_image()')
    return preview_path,preview_filename,keep_going


def open_image(imagename):
    hdu = fits.open(imagename)
    hdu = montage.reproject_hdu(hdu[0], north_aligned=True) 
    image = hdu.data
    nans = np.isnan(image)
    image[nans] = 0
    header = hdu.header
    wcs = WCS(header)
    return image, header, wcs


def main():

    # Define filenames and paths
    if len(sys.argv)!=5:
        print('Usage:\n python3 pyt_ingest_exposure_data_sdss.py [base_path] [sqlite_file] [thread_idx] [start_detection_id]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    thread_idx          = int(sys.argv[3])
    start_detection_id  = int(sys.argv[4])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('LN01_link_object_photometry_{:02d} execution started'.format(thread_idx),'LN01_link_object_photometry_{:02d} execution started.'.format(thread_idx))
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'LN01_link_object_photometry_{:02d}'.format(thread_idx))

        # Create and open detection data output files
        dir_detection_data  = base_path + 'detection_data/'
        preview_dir_path    = base_path + 'detection_previews/'
        preview_path_stem   = base_path + 'detection_previews/previews'
        preview_size_arcsec = 30
            
        if not os.path.isdir(dir_detection_data):
            mcd.create_directory(dir_detection_data,path_logfile,path_errorfile)
        if not os.path.isdir(preview_dir_path):
            mcd.create_directory(preview_dir_path,path_logfile,path_errorfile)
        output_data_file = dir_detection_data + 'detection_data_{:02d}_{:s}_toingest.txt'.format(thread_idx,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        with open(output_data_file,'w') as of:
            of.write('detection_id  x_coord_pred  y_coord_pred    x_coord    y_coord  right_ascension    declination  dist_pred_postn        flux_s       dflux_s        flux_t       dflux_t     mag_s  magerr_s     mag_t  magerr_t    fwhm_s        SNR  trail_len  trail_fwhm  trail_phi  flux_t_05rK  dflux_t_05rK  mag_t_05rK  magerr_t_05rK  flux_t_10rK  dflux_t_10rK  mag_t_10rK  magerr_t_10rK  flux_t_20rK  dflux_t_20rK  mag_t_20rK  magerr_t_20rK  flux_t_30rK  dflux_t_30rK  mag_t_30rK  magerr_t_30rK  flux_t_40rK  dflux_t_40rK  mag_t_40rK  magerr_t_40rK  flux_t_50rK  dflux_t_50rK  mag_t_50rK  magerr_t_50rK    ra_source_1   dec_source_1  dist_src_pred_1  dist_src_actual_1  mag_src_1  magerr_src_1    ra_source_2   dec_source_2  dist_src_pred_2  dist_src_actual_2  mag_src_2  magerr_src_2       ra_src_3      dec_src_3  dist_src_pred_3  dist_src_actual_3  mag_src_3  magerr_src_3  src_density_1arcmin  dist_edge_left  dist_edge_right  dist_edge_bottom  dist_edge_top  preview_image_file         preview_image_path\n')

        # Compute detection details
        link_object_photometry_multithread(thread_idx,start_detection_id,sqlite_file,output_data_file,preview_path_stem,preview_size_arcsec,path_logfile,path_errorfile)
        mcd.compress_file_gzip(output_data_file)
            
    mcd.send_status_email('LN01_link_object_photometry_{:02d} execution completed'.format(thread_idx),'LN01_link_object_photometry_{:02d} execution complete.'.format(thread_idx))
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()
