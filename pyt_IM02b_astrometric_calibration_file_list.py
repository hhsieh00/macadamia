### EDIT HISTORY ###
# 2021-06-19: combined SDSS and MegaCam code to make generic code module and added new logging functionality;
#              also now checks for interrupted processing and resets files for re-processing

import os, sys
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
import scipy
from scipy.stats import norm
import glob, bz2, subprocess
import os.path
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.wcs import WCS
from astropy.time import Time, TimeDelta
from astropy.io.votable import parse_single_table, VOTableSpecWarning, VOWarning
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from decimal import *
import sqlite3
from sqlite3 import Error
import jdcal
from jdcal import gcal2jd,jd2gcal
import statistics
from uncertainties import unumpy,ufloat
from uncertainties.umath import *
from uncertainties.unumpy import *
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE -- RETRIEVE DATA #####

def retrieve_instrument_id(sqlite_file,instrument_name,keep_going,path_logfile,path_errorfile):
    # Retrieve instrument_id corresponding to a given instrument name
    instrument_id = -1
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving instrument ID for {:s} from database...'.format(instrument_name))
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            query = "SELECT instrument_id FROM instruments WHERE instrument_name='{:s}'".format(instrument_name)
            mcd.output_log_entry(path_logfile,query)
            cursor.execute(query)
            row = cursor.fetchone()
            if row != None:
                instrument_id = int(row[0])
                mcd.output_log_entry(path_logfile,'instrument_id={:d}'.format(instrument_id))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,"Instrument '{:s}' not found.".format(instrument_name))
                mcd.send_status_email('IM02b - Function failed: retrieve_instrument_id()','IM02b - Function failed: retrieve_instrument_id()')
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: retrieve_instrument_id()')
            print(e)
            mcd.send_status_email('IM02b - Function failed: retrieve_instrument_id()','IM02b - Function failed: retrieve_instrument_id()')
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_instrument_id()')
    return(instrument_id,keep_going)


def retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve mosaic_element_id corresponding to a given mosaic element for a given instrument
    mosaic_element_id = -1
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving mosaic element ID from database...')
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            query = "SELECT mosaic_element_id FROM mosaic_elements WHERE instrument_id={:d} AND mosaic_element_num={:d}".format(instrument_id,mosaic_element_num)
            mcd.output_log_entry(path_logfile,query)
            cursor.execute(query)
            row = cursor.fetchone()
            if row != None:
                mosaic_element_id = int(row[0])
                mcd.output_log_entry(path_logfile,'Mosaic element {:d} (mosaic_element_id={:d}) for {:s} instrument found.'.format(mosaic_element_num,mosaic_element_id,instrument_name))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Mosaic element {:d} for {:s} instrument not found.'.format(mosaic_element_num,instrument_name))
                mcd.send_status_email('IM02b_{:d} - Function failed: retrieve_mosaic_element_id()'.format(thread_idx),'IM02b_{:d} - Function failed: retrieve_mosaic_element_id()'.format(thread_idx))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: retrieve_mosaic_element_id()')
            print(e)
            mcd.send_status_email('IM02b_{:d} - Function failed: retrieve_mosaic_element_id()'.format(thread_idx),'IM02b_{:d} - Function failed: retrieve_mosaic_element_id()'.format(thread_idx))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_mosaic_element_id()')
    return(mosaic_element_id,keep_going)


##### FUNCTION DEFINITIONS -- EXTRACT IMAGE METADATA #####

def process_date_tai_sdss(date_tai):
    date_tai_formatted = date_tai
    if date_tai[0:2].isnumeric() and date_tai[2] == '/' and date_tai[3:5].isnumeric() and date_tai[5] == '/' and date_tai[6:8].isnumeric:
        if int(date_tai[6:8]) > 20: year = '19' + date_tai[6:8]
        else: year = '20' + date_tai[6:8]
        date_tai_formatted = year + '-' + date_tai[3:5] + '-' + date_tai[0:2]
    return date_tai_formatted

def extract_header_info_sdss(fits_filename,processed_files_toingest_filename,keep_going,path_logfile,path_errorfile):
    # Extract header information from SDSS image file
    headerinfo = ['',0,'',0,'',0,0,0,0,0,1,1,0,0,0]
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Extracting header info from {:s}...'.format(fits_filename))
            hdulist = fits.open(fits_filename)
            hdr  = hdulist[0].header
            data = hdulist[0].data
            date_tai            = hdr['DATE-OBS']
            exposure_end_tai    = hdr['TAIHMS']
            exposure_time       = float(hdr['EXPTIME'])
            pointing_center_ra  = float(hdr['RA'])
            pointing_center_dec = float(hdr['DEC'])
            date_tai = process_date_tai_sdss(date_tai)
            t_end = Time('{:s} {:s}'.format(date_tai,exposure_end_tai),scale='tai',format='iso')
            dt = TimeDelta(exposure_time,format='sec')
            t_start = t_end - dt
            exposure_start_tai = t_start.iso[11:]
            exposure_start_jd  = float(t_start.jd)
            filter_name = fits_filename[6:7]
            airmass  = 1 / np.sin(float(hdr['ALT'])/180*math.pi)
            tracking_rate_ra  = float(0)
            tracking_rate_dec = float(0)
            binning_x = int(hdr['COLBIN'])
            binning_y = int(hdr['ROWBIN'])
            exp_npix_y,exp_npix_x = data.shape
            pixscale = 0.396
            mcd.output_log_entry(path_logfile,'Extracting header info done.')
            headerinfo = [date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra, \
                pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,pixscale]
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: extract_header_info_sdss()'.format(fits_filename))
            print(e)
            with open(fits_filename[:-5]+'.header_extraction_failed','w') as f:
                f.write('{:s} - Header info extraction using extract_header_info_sdss() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM02b - Function failed: extract_header_info_sdss()','IM02b - Function failed: extract_header_info_sdss()')
            with open(processed_files_toingest_filename,'a') as of:
                of.write('Header info extraction failed\n')
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,path_errorfile,'{:s} - Function skipped: extract_header_info_sdss()'.format(fits_filename))
    return headerinfo,keep_going


def extract_header_info_megacam(fits_filename,processed_files_toingest_filename,keep_going,path_logfile,path_errorfile):
    # Extract header information from MegaCam image file
    headerinfo = ['',0,'',0,'',0,0,0,0,0,1,1,0,0,0]
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Extracting header info from {:s}...'.format(fits_filename))
            hdulist = fits.open(fits_filename)
            hdr  = hdulist[0].header
            data = hdulist[0].data
            date_utc            = hdr['DATE-OBS']        # MegaCam FITS header comment: Date at start of observation (UTC)
            date_tai            = date_utc
            exposure_start_utc  = hdr['UTC-OBS']         # MegaCam FITS header comment: Time at start of observation (UTC)
            exposure_time       = float(hdr['EXPTIME'])  # MegaCam FITS header comment: Measured integration time (seconds)
            pointing_center_ra  = float(hdr['RA_DEG'])   # MegaCam FITS header comment: Object right ascension in degrees
            pointing_center_dec = float(hdr['DEC_DEG'])  # MegaCam FITS header comment: Object declination in degrees
            t_start = Time('{:s} {:s}'.format(date_utc,exposure_start_utc),scale='utc',format='iso')
            exposure_start_jd  = float(t_start.jd)
            exposure_start_tai_time = t_start.tai
            exposure_start_tai_time.format = 'iso'
            exposure_start_tai = exposure_start_tai_time.value[11:]
            pixscale = float(hdr['PIXSCAL1'])  # Pixel scale for axis 1 (arcsec/pixel)
            filter_id = hdr['FILTER']          # MegaCam FITS header comment: Name of filter in position
            if filter_id[0] == 'G' or filter_id[0] == 'g': filter_name = 'g'
            if filter_id[0] == 'R' or filter_id[0] == 'r': filter_name = 'r'
            if filter_id[0] == 'I' or filter_id[0] == 'i': filter_name = 'i'
            if filter_id[0] == 'Z' or filter_id[0] == 'z': filter_name = 'z'
            airmass  = float(hdr['AIRMASS'])   # MegaCam FITS header comment: Airmass at start of observation
            tracking_rate_ra  = float(0)
            tracking_rate_dec = float(0)
            binning_x = int(hdr['CCDBIN1'])    # MegaCam FITS header comment: Binning factor along first axis
            binning_y = int(hdr['CCDBIN2'])    # MegaCam FITS header comment: Binning factor along second axis
            exp_npix_y,exp_npix_x = data.shape
            mcd.output_log_entry(path_logfile,'Extracting header info done.')
            headerinfo = [date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra, \
                pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,pixscale]
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: extract_header_info_megacam()'.format(fits_filename))
            print(e)
            with open(fits_filename[:-5]+'.header_extraction_failed','w') as f:
                f.write('{:s} - Header info extraction using extract_header_info_megacam() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM02b - Function failed: extract_header_info_megacam()','IM02b - Function failed: extract_header_info_megacam()')
            with open(processed_files_toingest_filename,'a') as of:
                of.write('Header info extraction failed\n')
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: extract_header_info_megacam()')
    return headerinfo,keep_going


##### FUNCTION DEFINITIONS -- ASTROMETRIC CALIBRATION #####

def compute_astrometric_solution(fits_filename,processed_files_toingest_filename,pointing_ctr_ra,pointing_ctr_dec,pixscale,keep_going,path_logfile,path_errorfile):
    wcs_results = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Computing astrometric solution for {:s}...'.format(fits_filename))
            # Run astrometry.net code
            str_filename = fits_filename[:-5]
            cmd_solve = 'solve-field'
            if os.path.exists('/home/hhsieh/astrometry_gaia_pohaku.cfg'):
                cmd_config_file,value_config_file   = '--config','/home/hhsieh/astrometry_gaia_pohaku.cfg'
            elif os.path.exists('/home/hhsieh/astrometry_gaia_atlas.cfg'):
                cmd_config_file,value_config_file   = '--config','/home/hhsieh/astrometry_gaia_atlas.cfg'
            elif os.path.exists('/usr/local/astrometry/etc/astrometry.cfg'):
                cmd_config_file,value_config_file   = '--config','/usr/local/astrometry/etc/astrometry.cfg'
            else:
                cmd_config_file,value_config_file   = '--config',''
            cmd_ra_prediction,value_ra_prediction   = '--ra','{:f}'.format(pointing_ctr_ra)
            cmd_dec_prediction,value_dec_prediction = '--dec','{:f}'.format(pointing_ctr_dec)
            cmd_search_radius,value_search_radius   = '--radius','2'
            cmd_scale_units,value_scale_units       = '--scale-units','arcsecperpix'
            cmd_scale_low,value_scale_low           = '--scale-low','{:f}'.format(pixscale-0.1)
            cmd_scale_high,value_scale_high         = '--scale-high','{:f}'.format(pixscale+0.1)
            cmd_cpu_limit,value_cpu_limit           = '--cpulimit','60'
            if value_config_file != '':
                cmd = [cmd_solve,cmd_config_file,value_config_file,cmd_ra_prediction,value_ra_prediction,cmd_dec_prediction,value_dec_prediction,cmd_search_radius,value_search_radius,cmd_scale_units,value_scale_units,cmd_scale_low,value_scale_low,cmd_scale_high,value_scale_high,cmd_cpu_limit,value_cpu_limit,fits_filename]
                cmd_str = '{:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(cmd_solve,cmd_config_file,value_config_file,cmd_ra_prediction,value_ra_prediction,cmd_dec_prediction,value_dec_prediction,cmd_search_radius,value_search_radius,cmd_scale_units,value_scale_units,cmd_scale_low,value_scale_low,cmd_scale_high,value_scale_high,cmd_cpu_limit,value_cpu_limit,fits_filename)
            else:
                cmd = [cmd_solve,cmd_ra_prediction,value_ra_prediction,cmd_dec_prediction,value_dec_prediction,cmd_search_radius,value_search_radius,cmd_scale_units,value_scale_units,cmd_scale_low,value_scale_low,cmd_scale_high,value_scale_high,cmd_cpu_limit,value_cpu_limit,fits_filename]
                cmd_str = '{:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(cmd_solve,cmd_ra_prediction,value_ra_prediction,cmd_dec_prediction,value_dec_prediction,cmd_search_radius,value_search_radius,cmd_scale_units,value_scale_units,cmd_scale_low,value_scale_low,cmd_scale_high,value_scale_high,cmd_cpu_limit,value_cpu_limit,fits_filename)
            mcd.output_log_entry(path_logfile,cmd_str)
            process = subprocess.call(cmd)
            if os.path.isfile(str_filename + '.solved'):
                # Rename WCS-corrected FITS file
                if os.path.exists(str_filename+'_wcs.fits'):
                    datetime_str = datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
                    os.rename(str_filename+'_wcs.fits',str_filename+'_wcs.{:s}.fits'.format(datetime_str))
                    mcd.compress_file_fpack(str_filename+'_wcs.{:s}.fits'.format(datetime_str))
                os.rename(str_filename+'.new',str_filename+'_wcs.fits')
                # Extract number of sources used to derive WCS solution
                filename_rdls = str_filename + '.rdls'
                hdulist = fits.open(filename_rdls)
                wcs_nstars = len(hdulist[1].data)
                hdulist.close()
                # Extract WCS results
                hdulist = fits.open(str_filename+'_wcs.fits')
                hdr = hdulist[0].header
                wcs_crpix1 = float(hdr['CRPIX1'])
                wcs_crpix2 = float(hdr['CRPIX2'])
                wcs_crval1 = float(hdr['CRVAL1'])
                wcs_crval2 = float(hdr['CRVAL2'])
                wcs_cd1_1  = float(hdr['CD1_1'])
                wcs_cd1_2  = float(hdr['CD1_2'])
                wcs_cd2_1  = float(hdr['CD2_1'])
                wcs_cd2_2  = float(hdr['CD2_2'])
                hdulist.close()
                wcs_pixscale_x,wcs_pixscale_y = compute_wcs_pixel_scale(str_filename+'_wcs.fits',path_logfile,path_errorfile)
                wcs_results = [wcs_nstars,wcs_crpix1,wcs_crpix2,wcs_crval1,wcs_crval2,wcs_cd1_1,wcs_cd1_2,wcs_cd2_1,wcs_cd2_2,wcs_pixscale_x,wcs_pixscale_y]
                keep_going = True
                mcd.output_log_entry(path_logfile,'Astrometric solution successful.')
                with open(fits_filename[:-5]+'.wcssolved','w') as f:
                    f.write('{:s} - Astrometric calibration successful.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                with open(processed_files_toingest_filename,'a') as of:
                    of.write('Astrometric calibration successful\n')
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Astrometric solution unsuccessful.'.format(fits_filename))
                with open(fits_filename[:-5]+'.astrometric_calibration_failed','w') as f:
                    f.write('{:s} - Astrometric calibration using compute_astrometric_solution() unsuccessful.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
                with open(processed_files_toingest_filename,'a') as of:
                    of.write('Astrometric calibration unsuccessful\n')
                mcd.compress_file_fpack(fits_filename)
                os.remove(fits_filename[:-5]+'-objs.png')
                os.remove(fits_filename[:-5]+'.axy')
                mcd.send_status_email('IM02b - Astrometric calibration failed for {:s}'.format(fits_filename),'Astrometric calibration for {:s} using compute_astrometric_solution() unsuccessful.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),fits_filename))
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: compute_astrometric_solution()'.format(fits_filename))
            print(e)
            with open(fits_filename[:-5]+'.astrometric_calibration_failed','w') as f:
                f.write('{:s} - Astrometric calibration using compute_astrometric_solution() failed.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            with open(processed_files_toingest_filename,'a') as of:
                of.write('Astrometric calibration unsuccessful\n')
            mcd.send_status_email('IM02b - Function failed: compute_astrometric_solution()','IM02b - Function failed: compute_astrometric_solution()'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_astrometric_solution()')
    return wcs_results,keep_going


def compute_wcs_pixel_scale(fits_file,path_logfile,path_errorfile):
    wcs_pixscale_x,wcs_pixscale_y = -1,-1
    try:
        ra1,dec1 = mcd.get_radec_from_pixel_coords(fits_file,1,1)
        ra2,dec2 = mcd.get_radec_from_pixel_coords(fits_file,1001,1)
        ra3,dec3 = mcd.get_radec_from_pixel_coords(fits_file,1,1001)
        wcs_pixscale_x  = mcd.compute_angular_dist_radec(ra1,dec1,ra2,dec2)/1000*3600
        wcs_pixscale_y  = mcd.compute_angular_dist_radec(ra1,dec1,ra3,dec3)/1000*3600
    except:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: compute_wcs_pixel_scale()'.format(fits_file))
    return wcs_pixscale_x,wcs_pixscale_y


##### FUNCTION DEFINITIONS -- DATA EXPORTING #####

def write_astrometry_data_toingest(instrument,base_path,dir_raw_data,dir_processed_data,sub_dir,data_dir,astrometry_data_toingest_filename,fits_filestem,mosaic_elem_table,headerinfo,wcs_results,keep_going,path_logfile,path_errorfile):
    # writes table of exposure data for later ingestion into database
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Writing exposure data to {:s} for ingestion...'.format(astrometry_data_toingest_filename))
            with open(astrometry_data_toingest_filename,'a') as f:
                filename_base = fits_filestem[:-4]       # remove extension number from FITS filename
                ext_number    = int(fits_filestem[-3:])  # get extension number from last 3 digits of fits_filestem
                mosaic_element_id = -1
                for idx in range(0,len(mosaic_elem_table)):
                    if mosaic_elem_table[idx]['elem_idx'] == ext_number:
                        mosaic_element_id = mosaic_elem_table[idx]['element_id']                        
                if instrument == 'SDSS':
                    filename_parts      = filename_base.split('-')
                    raw_data_link       = 'https://data.sdss.org/sas/dr9/boss/photoObj/frames/301/{:d}/{:d}/{:s}.fits.bz2'.format(int(filename_parts[2]),int(filename_parts[3]),filename_base)
                    raw_data_path       = dir_raw_data+sub_dir+data_dir
                    proc_data_path      = dir_processed_data+sub_dir+data_dir
                if instrument == 'MegaCam':
                    raw_data_link       = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/CFHT/{:s}.fits.fz'.format(filename_base)
                    raw_data_path       = dir_raw_data+data_dir
                    proc_data_path      = dir_processed_data+data_dir
                base_filename       = filename_base+'.fits'
                raw_data_file       = filename_base+'.fits.fz'
                proc_data_file      = fits_filestem+'_wcs.fits.fz'
                date_tai            = headerinfo[0]
                exposure_start_jd   = float(headerinfo[1])
                exposure_start_tai  = headerinfo[2]
                exposure_time       = float(headerinfo[3])
                filter_name         = headerinfo[4]
                airmass             = float(headerinfo[5])
                pointing_center_ra  = float(headerinfo[6])
                pointing_center_dec = float(headerinfo[7])
                tracking_rate_ra    = float(headerinfo[8])
                tracking_rate_dec   = float(headerinfo[9])
                binning_x           = int(headerinfo[10])
                binning_y           = int(headerinfo[11])
                exp_npix_x          = int(headerinfo[12])
                exp_npix_y          = int(headerinfo[13])
                wcs_nstars          = int(wcs_results[0])
                wcs_crpix1          = int(wcs_results[1])
                wcs_crpix2          = int(wcs_results[2])
                wcs_crval1          = float(wcs_results[3])
                wcs_crval2          = float(wcs_results[4])
                wcs_cd1_1           = float(wcs_results[5])
                wcs_cd1_2           = float(wcs_results[6])
                wcs_cd2_1           = float(wcs_results[7])
                wcs_cd2_2           = float(wcs_results[8])
                wcs_pixscale_x      = float(wcs_results[9])
                wcs_pixscale_y      = float(wcs_results[10])
                f.write('{:>17d}   {:<27s}   {:<10s}   {:18.10f}         {:s}  {:15.10f}  {:<11s}   {:7.3f}   {:15.10f}   {:16.10f}   {:16.10f}   {:17.10f}   {:>5d}   {:>5d}   {:>10d}   {:>10d}   {:>10d}   {:>10d}   {:>10d}   {:15.10f}   {:15.10f}   {:13.10f}   {:13.10f}   {:13.10f}   {:13.10f}   {:14.10f}   {:14.10f}   {:s}   {:s}   {:s}   {:s}   {:s}\n'.format(mosaic_element_id,base_filename,date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra,pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,wcs_nstars,wcs_crpix1,wcs_crpix2,wcs_crval1,wcs_crval2,wcs_cd1_1,wcs_cd1_2,wcs_cd2_1,wcs_cd2_2,wcs_pixscale_x,wcs_pixscale_y,raw_data_link,raw_data_path,raw_data_file,proc_data_path,proc_data_file))
            mcd.output_log_entry(path_logfile,'Writing exposure data done.')
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: write_astrometry_data_toingest()'.format(fits_filestem+'.fits'))
            print(e)
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: write_astrometry_data_toingest()')
    return keep_going


##### FUNCTION DEFINITIONS -- CLEAN UP FILES #####

def clean_previous_wcs_output_file(prev_wcs_output_file):
    # Deletes specified astrometry.net output file or gzip-compressed equivalent, if present
    if os.path.exists(prev_wcs_output_file):
        os.remove(prev_wcs_output_file)
    if os.path.exists(prev_wcs_output_file+'.gz'):
        os.remove(prev_wcs_output_file+'.gz')
    return None

def clean_previous_wcs_output_files(fits_filename,orig_data_path):
    # Removes all standard astrometry.net output
    clean_previous_wcs_output_file(fits_filename[:-5]+'.axy')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.corr')
    clean_previous_wcs_output_file(fits_filename[:-5]+'-indx.png')
    clean_previous_wcs_output_file(fits_filename[:-5]+'-indx.xyls')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.match')
    clean_previous_wcs_output_file(fits_filename[:-5]+'-ngc.png')
    clean_previous_wcs_output_file(fits_filename[:-5]+'-objs.png')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.rdls')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.solved')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.wcs')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.new')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.wcssolved')
    clean_previous_wcs_output_file(fits_filename[:-5]+'.astrometric_calibration_failed')
    clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs_wcs.fits')
    clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs_wcs.fits.fz')
    clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.wcssolved')
    clean_previous_wcs_output_file(orig_data_path+fits_filename[:-5]+'_wcs.fits.fz')
    return None

def check_clean_previous_wcs_output(data_path,fits_filename,keep_going,path_logfile,path_errorfile):
    # Checks for presence of previous astrometry.net output indicating previously interrupted processing
    # If previous astrometry.net output found, delete previous output and reset FITS file for re-processing
    if keep_going:
        wcs_output_path = data_path+'wcs_output/'
        orig_data_path  = data_path+'original_data/'
        mcd.create_directory(wcs_output_path,path_logfile,path_errorfile)
        mcd.create_directory(orig_data_path,path_logfile,path_errorfile)
        if os.path.exists(fits_filename[:-5]+'.axy')       or os.path.exists(fits_filename[:-5]+'.axy.gz') \
        or os.path.exists(fits_filename[:-5]+'.corr')      or os.path.exists(fits_filename[:-5]+'.corr.gz') \
        or os.path.exists(fits_filename[:-5]+'-indx.png')  or os.path.exists(fits_filename[:-5]+'-indx.png.gz') \
        or os.path.exists(fits_filename[:-5]+'-indx.xyls') or os.path.exists(fits_filename[:-5]+'-indx.xyls.gz') \
        or os.path.exists(fits_filename[:-5]+'.match')     or os.path.exists(fits_filename[:-5]+'.match.gz') \
        or os.path.exists(fits_filename[:-5]+'-ngc.png')   or os.path.exists(fits_filename[:-5]+'-ngc.png.gz') \
        or os.path.exists(fits_filename[:-5]+'-objs.png')  or os.path.exists(fits_filename[:-5]+'-objs.png.gz') \
        or os.path.exists(fits_filename[:-5]+'.rdls')      or os.path.exists(fits_filename[:-5]+'.rdls.gz') \
        or os.path.exists(fits_filename[:-5]+'.solved')    or os.path.exists(fits_filename[:-5]+'.solved.gz') \
        or os.path.exists(fits_filename[:-5]+'.wcs')       or os.path.exists(fits_filename[:-5]+'.wcs.gz') \
        or os.path.exists(fits_filename[:-5]+'.new')       or os.path.exists(fits_filename[:-5]+'_wcs.fits'):
            # Previous intermediate astrometry.net output found
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Previous intermediate astrometry.net output found for {:s}...'.format(fits_filename))
            if os.path.exists(fits_filename[:-5]+'.fits.fz'):
                # Original compressed FITS file found in main directory
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} found; removing other intermediate astrometry.net output...'.format(fits_filename[:-5]+'.fits.fz'))
                clean_previous_wcs_output_files(fits_filename,orig_data_path)
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits')
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits.fz')
            elif os.path.exists(fits_filename[:-5]+'.fits'):
                # Original uncompressed FITS file found in main directory
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} found; compressing file and removing other intermediate astrometry.net output...'.format(fits_filename[:-5]+'.fits'))
                mcd.compress_file_fpack(fits_filename[:-5]+'.fits')
                clean_previous_wcs_output_files(fits_filename,orig_data_path)
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits')
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits.fz')
            elif os.path.exists(orig_data_path+fits_filename[:-5]+'.fits.fz'):
                # Original compressed FITS file found in original data directory
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} found; returning to main directory and removing other intermediate astrometry.net output...'.format(orig_data_path+fits_filename[:-5]+'.fits.fz'))
                mcd.compress_file_fpack(fits_filename[:-5]+'.fits')
                clean_previous_wcs_output_files(fits_filename,orig_data_path)
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits')
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits.fz')
            elif os.path.exists(fits_filename[:-5]+'_wcs.fits.fz'):
                # Original FITS file not found, but compressed WCS-solved FITS file found
                # Rename to original FITS filename for reprocessing
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} found but {:s} not found'.format(fits_filename[:-5]+'_wcs.fits.fz',fits_filename[:-5]+'.fits.fz'))
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Renaming {:s} to {:s}; removing other intermediate astrometry.net output...'.format(fits_filename[:-5]+'_wcs.fits.fz',fits_filename[:-5]+'.fits.fz'))
                os.rename(fits_filename[:-5]+'_wcs.fits.fz',fits_filename[:-5]+'.fits.fz')
                clean_previous_wcs_output_files(fits_filename,orig_data_path)
                clean_previous_wcs_output_file(fits_filename[:-5]+'_wcs.fits')
            elif os.path.exists(fits_filename[:-5]+'_wcs.fits'):
                # Original FITS file not found, but uncompressed WCS-solved FITS file found
                # Rename to original FITS filename for reprocessing
                mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} found but {:s} not found'.format(fits_filename[:-5]+'_wcs.fits',fits_filename[:-5]+'.fits'))
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Renaming {:s} to {:s}; compressing file and removing other intermediate astrometry.net output...'.format(fits_filename[:-5]+'_wcs.fits',fits_filename[:-5]+'.fits'))
                os.rename(fits_filename[:-5]+'_wcs.fits',fits_filename[:-5]+'.fits')
                mcd.compress_file_fpack(fits_filename[:-5]+'.fits')
                clean_previous_wcs_output_files(fits_filename,orig_data_path)
    return keep_going

def clean_wcs_output_file(wcs_output_file,wcs_output_path):
    # Checks wcs output directory for previous astrometry.net output; removes previous output if found
    # Then compresses current astrometry.net output file and moves to wcs output directory
    if os.path.exists(wcs_output_file):
        if os.path.exists(wcs_output_file+'.gz'):                 os.remove(wcs_output_file+'.gz')
        if os.path.exists(wcs_output_path+wcs_output_file+'.gz'): os.remove(wcs_output_path+wcs_output_file+'.gz')
        mcd.compress_file_gzip(wcs_output_file)
        os.rename(wcs_output_file+'.gz',wcs_output_path+wcs_output_file+'.gz')
    return None

def clean_wcs_output(path_datadownload,fits_filename,keep_going,path_logfile,path_errorfile):
    # Archives wcs output in wcs output directory
    if keep_going:
        mcd.output_log_entry(path_logfile,'Cleaning up astrometry.net output...')
        wcs_output_path = path_datadownload+'wcs_output/'
        orig_data_path = path_datadownload+'original_data/'
        mcd.create_directory(wcs_output_path,path_logfile,path_errorfile)
        mcd.create_directory(orig_data_path,path_logfile,path_errorfile)
        try:
            clean_wcs_output_file(fits_filename[:-5]+'.axy',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'.corr',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'-indx.png',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'-indx.xyls',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'.match',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'-ngc.png',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'-objs.png',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'.rdls',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'.solved',wcs_output_path)
            clean_wcs_output_file(fits_filename[:-5]+'.wcs',wcs_output_path)
            if os.path.exists(fits_filename):
                if os.path.exists(fits_filename+'.fz'):                os.remove(fits_filename+'.fz')
                if os.path.exists(orig_data_path+fits_filename+'.fz'): os.remove(orig_data_path+fits_filename+'.fz')
                mcd.compress_file_fpack(fits_filename)
                os.rename(fits_filename+'.fz',orig_data_path+fits_filename+'.fz')
            mcd.output_log_entry(path_logfile,'Cleaning up astrometry.net output done.')
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: clean_wcs_output()')
            print(e)
            with open(fits_filename[:-5]+'.wcs_output_cleanup_failed','w') as f:
                f.write('{:s} - WCS output cleaning using clean_wcs_output() failed.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM02b - Function failed: clean_wcs_output()','IM02b - Function failed: clean_wcs_output()')
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: clean_wcs_output()')
    return keep_going


def main():

    # Define filenames and paths
    if len(sys.argv)!=5:
        print('Usage:\n python3 pyt_IM02b_astrometric_calibration_file_list.py [base_path] [sqlite_file] [instrument] [thread_idx]\n')
        print(' (Trailing /\ needed in path specification)\n')
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    instrument  = sys.argv[3]
    thread_idx  = int(sys.argv[4])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('IM02_astrometric_calibration_file_list_{:s}_{:d} execution started.'.format(instrument,thread_idx),'IM02_astrometric_calibration_file_list_{:s}_{:02d} execution started.'.format(instrument,thread_idx))

        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'IM02b_astrometric_calibration_file_list_{:s}_{:02d}'.format(instrument,thread_idx))
        
        os.chdir(base_path)
        if instrument == 'SDSS':
            dir_raw_data       = base_path+'data_sdss_raw/'
            dir_processed_data = base_path+'data_sdss_processed/'
        elif instrument == 'MegaCam':
            dir_raw_data       = base_path+'data_megacam_raw/'
            dir_processed_data = base_path+'data_megacam_processed/'
        else:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Instrument {:s} not recognized or not yet implemented.'.format(instrument))
            keep_going = False
            
    if keep_going:
        dir_exposure_data  = base_path+'exposure_data/'
        dir_files_toproc   = base_path+'ssois_queries/files_toproc/'
        mcd.create_directory(dir_exposure_data,path_logfile,path_errorfile)
        if not os.path.isdir(dir_raw_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Raw data directory {:s} not found.'.format(dir_raw_data))
            mcd.send_status_email('IM02_astrometric_calibration_sdss_file_list_{:s}_{:02d} execution failed'.format(instrument,thread_idx),'IM02b_astrometric_photometric_calibration_sdss_file_list_{:s}_{:02d} execution failed - Raw data directory {:s} not found.'.format(instrument,thread_idx,dir_raw_data))
            keep_going = False
        if not os.path.isdir(dir_processed_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Processed data directory {:s} not found.'.format(dir_processed_data))
            mcd.send_status_email('IM02_astrometric_calibration_sdss_file_list_{:s}_{:02d} execution failed'.format(instrument,thread_idx),'IM02b_astrometric_photometric_calibration_sdss_file_list_{:s}_{:02d} execution failed - Processed data directory {:s} not found.'.format(instrument,thread_idx,dir_processed_data))
            keep_going = False
        if not os.path.isdir(dir_files_toproc):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'File lists directory {:s} not found.'.format(dir_files_toproc))
            mcd.send_status_email('IM02_astrometric_calibration_sdss_file_list_{:s}_{:02d} execution failed'.format(instrument,thread_idx),'IM02b_astrometric_photometric_calibration_sdss_file_list_{:s}_{:02d} execution failed - File lists directory {:s} not found.'.format(instrument,thread_idx,dir_files_toproc))
            keep_going = False
    
    if keep_going:
        match_threshold_pix = 3
        
        if instrument == 'SDSS':
            # SDSS-specific data
            instrument_name = 'SDSS Imaging Camera'
            first_element = 0
            num_elements  = 1
            mosaic_elem_datarows = [(0,-1)]  # element index number and dummy value for element_id

            # Retrieve telescope/instrument data
            #instrument_id,keep_going,log_file = retrieve_instrument_id(sqlite_file,instrument_name,keep_going,log_file)
            instrument_id = 1  # SDSS Imaging Camera
            mosaic_elem_table    = Table(rows=mosaic_elem_datarows,names=('elem_idx','element_id'))
            for idx in range(0,num_elements):
                mosaic_element_num = mosaic_elem_table[idx]['elem_idx']
                #mosaic_element_id,keep_going,log_file = retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,keep_going,log_file)
                mosaic_element_id = 1  # SDSS Imaging Camera - Ext 0
                mosaic_elem_table[idx]['element_id'] = mosaic_element_id
                
        if instrument == 'MegaCam':
            ### MegaCam-specific data ###
            instrument_name = 'MegaCam'
            first_element = 1
            num_elements  = 40
            mosaic_elem_datarows = [[0 for idx2 in range(2)] for idx1 in range(num_elements)]
            for idx in range(0,num_elements):
                mosaic_elem_datarows[idx][0] = first_element + idx
                
            # Retrieve telescope/instrument data
            #instrument_id,keep_going,log_file = retrieve_instrument_id(sqlite_file,instrument_name,keep_going,log_file)
            instrument_id = 2  # MegaCam
            mosaic_elem_table = Table(rows=mosaic_elem_datarows,names=('elem_idx','element_id'))
            for idx in range(0,num_elements):
                mosaic_element_num = mosaic_elem_table[idx]['elem_idx']
                mosaic_element_id,keep_going = retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile)
                mosaic_elem_table[idx]['element_id'] = mosaic_element_id
                
    if keep_going:
        # Process data
        os.chdir(dir_processed_data)
        for ssois_query_results_filepath_gz in sorted(glob.glob(dir_files_toproc+'*_{:s}_{:02d}.txt.gz'.format(instrument,thread_idx))):
            mcd.output_log_entry(path_logfile,'Processing image files in {:s}...'.format(ssois_query_results_filepath_gz))
            # initialize exposure data output files
            astrometry_data_toingest_filename = dir_exposure_data+'00_astrometry_data_{:s}_{:s}_{:02d}_toingest.txt'.format(instrument,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),thread_idx)
            processed_files_toingest_filename = dir_exposure_data+'00_processed_files_astrometry_{:s}_{:s}_{:02d}_toingest.txt'.format(instrument,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),thread_idx)
            mcd.output_log_entry(path_logfile,'Writing astrometric output to {:s}...'.format(astrometry_data_toingest_filename))
            mcd.output_log_entry(path_logfile,'Writing processing statuses to {:s}...'.format(processed_files_toingest_filename))
            if not os.path.exists(astrometry_data_toingest_filename):
                with open(astrometry_data_toingest_filename,'w') as f:
                    f.write('mosaic_element_id   base_filename                   date_tai    exposure_start_jd   exposure_start_tai    exposure_time  filter_name   airmass   pointing_ctr_ra   pointing_ctr_dec   tracking_rate_ra   tracking_rate_dec   bin_x   bin_y   exp_npix_x   exp_npix_y   wcs_nstars   wcs_crpix1   wcs_crpix2        wcs_crval1        wcs_crval2       wcs_cd1_1       wcs_cd1_2       wcs_cd2_1       wcs_cd2_2   wcs_pixscale_x   wcs_pixscale_y   raw_data_link   raw_data_path   raw_data_file   proc_data_path   proc_data_file\n')
            if not os.path.exists(processed_files_toingest_filename):
                with open(processed_files_toingest_filename,'w') as f:
                    f.write('proc_data_file                                       status\n')
                
            # read in query results file data
            mcd.decompress_file_gzip(ssois_query_results_filepath_gz)
            ssois_query_results_filepath = ssois_query_results_filepath_gz[:-3]
            with open(ssois_query_results_filepath,'r') as ssois_query_results_file:
                for _ in range(1): #skip first 1 header line
                    next(ssois_query_results_file)
                for line in ssois_query_results_file:
                    image_link = line[136:-1].strip()
                    filename_toprocess,filename_found,idx = '',False,len(image_link)-1
                    while not filename_found:
                        if image_link[idx] != '/':
                            filename_toprocess = image_link[idx] + filename_toprocess
                        else:
                            filename_found = True
                        idx -= 1
                    if instrument == 'SDSS':
                        sub_dir  = '{:d}/'.format(int(filename_toprocess[8:14]))
                        data_dir = '{:d}/'.format(int(filename_toprocess[15]))
                        data_path = dir_processed_data+sub_dir+data_dir
                    if instrument == 'MegaCam':
                        obs_date   = line[12:22]
                        data_dir   = 'ut'+obs_date[:4]+obs_date[5:7]+obs_date[8:]
                        data_path = dir_processed_data+data_dir+'/'
                        
                    if os.path.isdir(data_path):
                        os.chdir(data_path)
                        proc_files_found,files_to_proc_found = False,False
                        for fits_filename_fz in sorted(glob.glob(filename_toprocess[:-9]+'_???.fits.fz')):
                            fits_filename = fits_filename_fz[:-3]  # remove .fz extension from filename
                            fits_filestem = fits_filename[:-5]     # remove .fits extension from filename
                            check_clean_previous_wcs_output(data_path,fits_filename,keep_going,path_logfile,path_errorfile)
                            if not os.path.exists(fits_filename_fz[:-8]+'.astrometric_calibration_failed') and not os.path.exists(fits_filename_fz[:-8]+'.wcssolved'):
                                keep_going = True
                                files_to_proc_found = True
                                mcd.output_log_entry(path_logfile,'Decompressing {:s}{:s} for processing...'.format(data_path,fits_filename_fz))
                                with open(processed_files_toingest_filename,'a') as of:
                                    of.write('{:<50s}   '.format(fits_filename_fz))
                                mcd.decompress_file_funpack(fits_filename_fz)
                
                                ### Extract header info from FITS file
                                if instrument == 'SDSS':
                                    headerinfo,keep_going = extract_header_info_sdss(fits_filename,processed_files_toingest_filename,keep_going,path_logfile,path_errorfile)
                                if instrument == 'MegaCam':
                                    headerinfo,keep_going = extract_header_info_megacam(fits_filename,processed_files_toingest_filename,keep_going,path_logfile,path_errorfile)
                                exptime          = float(headerinfo[3])
                                filter_name      = headerinfo[4]
                                pointing_ctr_ra  = headerinfo[6]
                                pointing_ctr_dec = headerinfo[7]
                                exp_npix_x       = headerinfo[12]
                                exp_npix_y       = headerinfo[13]
                                pixscale         = headerinfo[14]
                
                                ### Compute astrometric solution for FITS file
                                wcs_results,keep_going = compute_astrometric_solution(fits_filename,processed_files_toingest_filename,pointing_ctr_ra,pointing_ctr_dec,pixscale,keep_going,path_logfile,path_errorfile)
                                ra_deg     = wcs_results[3]
                                dec_deg    = wcs_results[4]
                                wcs_pixscale_x = wcs_results[9]
                                wcs_pixscale_y = wcs_results[10]
                                keep_going = clean_wcs_output(data_path,fits_filename,keep_going,path_logfile,path_errorfile)
                                fits_filename_wcs = fits_filestem + '_wcs.fits'  # e.g., file_basename_000 --> file_basename_000_wcs.fits
                                        
                                keep_going = write_astrometry_data_toingest(instrument,base_path,dir_raw_data,dir_processed_data,sub_dir,data_dir,astrometry_data_toingest_filename,fits_filestem,mosaic_elem_table,headerinfo,wcs_results,keep_going,path_logfile,path_errorfile)
                                
                                # Compress fits file for archiving
                                mcd.output_log_entry(path_logfile,'Compressing {:s}{:s}...'.format(data_path,fits_filename_wcs))
                                mcd.compress_file_fpack(fits_filename_wcs)
                        if not files_to_proc_found:
                            mcd.output_log_entry(path_logfile,'No extensions found requiring processing for {:s}{:s}.'.format(data_path,filename_toprocess))
                    else:
                        mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Data path {:s} not found.'.format(filename_toprocess,data_path))
                        mcd.send_status_email('IM02b_astrometric_calibration_file_list_{:s}_{:02d} execution failed'.format(instrument,thread_idx),'IM02b_astrometric_calibration_file_list_{:s}_{:02d} execution failed - {:s} - Data path {:s} not found.'.format(instrument,thread_idx,filename_toprocess,data_path))                            
                        
            mcd.remove_output_file_if_empty(path_logfile,astrometry_data_toingest_filename)
            mcd.remove_output_file_if_empty(path_logfile,processed_files_toingest_filename)
            
            if os.path.exists(astrometry_data_toingest_filename): mcd.compress_file_gzip(astrometry_data_toingest_filename)
            if os.path.exists(processed_files_toingest_filename): mcd.compress_file_gzip(processed_files_toingest_filename)
            mcd.compress_file_gzip(ssois_query_results_filepath)
                        
    mcd.send_status_email('IM02b_astrometric_calibration_file_list_{:s}_{:02d} execution complete.'.format(instrument,thread_idx),'IM02b_astrometric_calibration_file_list_{:s}_{:02d} execution complete.'.format(instrument,thread_idx))

    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

