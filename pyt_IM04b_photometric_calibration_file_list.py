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

### Edit History ###
# 2021-06-30: Added new logging functionality; made separate function for multiap photometry for each aperture
# 2021-06-30: Removed per-function notification emails for individual image processing failures


##### FUNCTION DEFINITIONS -- DATABASE -- RETRIEVE DATA #####

def retrieve_instrument_id(sqlite_file,instrument_name,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve instrument_id corresponding to a given instrument name
    instrument_id = -1
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving instrument ID for {:s} from database...'.format(instrument_name))
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            query = "SELECT instrument_id FROM instruments WHERE instrument_name='{:s}'".format(instrument_name)
            cursor.execute(query)
            mcd.output_log_entry(path_logfile,query)
            row = cursor.fetchone()
            if row != None:
                instrument_id = int(row[0])
                mcd.output_log_entry(path_logfile,'instrument_id={:d}'.format(instrument_id))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Instrument {:s} not found.'.format(instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: retrieve_instrument_id()'.format(instrument_name))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_instrument_id()')
    return instrument_id,keep_going


def retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve mosaic_element_id corresponding to a given mosaic element for a given instrument
    mosaic_element_id = -1
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving mosaic element ID from database...')
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            query = "SELECT mosaic_element_id FROM mosaic_elements WHERE instrument_id={:d} AND mosaic_element_num={:d}".format(instrument_id,mosaic_element_num)
            cursor.execute(query)
            mcd.output_log_entry(path_logfile,query)
            row = cursor.fetchone()
            if row != None:
                mosaic_element_id = int(row[0])
                mcd.output_log_entry(path_logfile,'Mosaic element {:d} (mosaic_element_id={:d}) found for {:s}.'.format(mosaic_element_num,mosaic_element_id,instrument_name))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Mosaic element {:d} for {:s} not found.'.format(mosaic_element_num,instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}, Ext {:d}: retrieve_mosaic_element_id()'.format(instrument_name,mosaic_element_num))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_mosaic_element_id()')
    return mosaic_element_id,keep_going


##### FUNCTION DEFINITIONS -- EXTRACT IMAGE METADATA #####

def process_date_tai_sdss(date_tai):
    date_tai_formatted = date_tai
    if date_tai[0:2].isnumeric() and date_tai[2] == '/' and date_tai[3:5].isnumeric() and date_tai[5] == '/' and date_tai[6:8].isnumeric:
        if int(date_tai[6:8]) > 20: year = '19' + date_tai[6:8]
        else: year = '20' + date_tai[6:8]
        date_tai_formatted = year + '-' + date_tai[3:5] + '-' + date_tai[0:2]
    return date_tai_formatted

def extract_header_info_sdss(fits_filename,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Extract header information from SDSS image file
    headerinfo = ['',0,'',0,'',0,0,0,0,0,1,1,0,0,0,0]
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
            ra_deg    = float(hdr['CRVAL1'])
            dec_deg   = float(hdr['CRVAL2'])
            exp_npix_y,exp_npix_x = data.shape
            mcd.output_log_entry(path_logfile,'Extracting header info done.')
            headerinfo = [date_tai,exposure_start_jd,exposure_start_tai,exposure_time,filter_name,airmass,pointing_center_ra, \
                pointing_center_dec,tracking_rate_ra,tracking_rate_dec,binning_x,binning_y,exp_npix_x,exp_npix_y,ra_deg,dec_deg]
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: extract_header_info_sdss()'.format(fits_filename))
            print(e)
            if os.path.exists(fits_filename[:-5]+'.photfailed'):
                with open(fits_filename[:-5]+'.photfailed','a') as f:
                    f.write('{:s} - SDSS header info extraction using extract_header_info_sdss() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename[:-5]+'.photfailed','w') as f:
                    f.write('{:s} - SDSS header info extraction using extract_header_info_sdss() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'SDSS header info extraction using extract_header_info_sdss() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: extract_header_info_sdss()')
    return headerinfo,proc_status,keep_going


##### FUNCTION DEFINITIONS -- ASTROMETRIC FUNCTIONS #####

def compute_wcs_pixel_scale(fits_file):
    wcs_pixscale_x,wcs_pixscale_y = -1,-1
    try:
        ra1,dec1 = mcd.get_radec_from_pixel_coords(fits_file,1,1)
        ra2,dec2 = mcd.get_radec_from_pixel_coords(fits_file,1001,1)
        ra3,dec3 = mcd.get_radec_from_pixel_coords(fits_file,1,1001)
        wcs_pixscale_x  = mcd.compute_angular_dist_radec(ra1,dec1,ra2,dec2)/1000*3600
        wcs_pixscale_y  = mcd.compute_angular_dist_radec(ra1,dec1,ra3,dec3)/1000*3600
    except:
        print('{:s} - Function failed: compute_wcs_pixel_scale()'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    return wcs_pixscale_x,wcs_pixscale_y


def get_radec_from_pixel_coords_array(fits_filename_wcs,field_source_table,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Convert pixel coordinates to RA/Dec using WCS information
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Converting pixel coordinates to RA/Dec using WCS information for {:s}...'.format(fits_filename_wcs))
            header = fits.getheader(fits_filename_wcs)
            w = WCS(header)
            ra_array,dec_array = w.wcs_pix2world(field_source_table['x_t'],field_source_table['y_t'],1)
            field_source_table['RA']          = ra_array
            field_source_table['Dec']         = dec_array
            field_source_table['ps1_RA']      = -999.0
            field_source_table['ps1_Dec']     = -999.0
            field_source_table['gp1_mag']     = -999.0
            field_source_table['gp1_err']     = -999.0
            field_source_table['rp1_mag']     = -999.0
            field_source_table['rp1_err']     = -999.0
            field_source_table['ip1_mag']     = -999.0
            field_source_table['ip1_err']     = -999.0
            field_source_table['zp1_mag']     = -999.0
            field_source_table['zp1_err']     = -999.0
            field_source_table['g_sdss_mag']  = -999.0
            field_source_table['g_sdss_err']  = -999.0
            field_source_table['r_sdss_mag']  = -999.0
            field_source_table['r_sdss_err']  = -999.0
            field_source_table['i_sdss_mag']  = -999.0
            field_source_table['i_sdss_err']  = -999.0
            field_source_table['z_sdss_mag']  = -999.0
            field_source_table['z_sdss_err']  = -999.0
            field_source_table['B_mag']       = -999.0
            field_source_table['B_err']       = -999.0
            field_source_table['V_mag']       = -999.0
            field_source_table['V_err']       = -999.0
            field_source_table['Rc_mag']      = -999.0
            field_source_table['Rc_err']      = -999.0
            field_source_table['Ic_mag']      = -999.0
            field_source_table['Ic_err']      = -999.0
            field_source_table['target_dist'] = -999.0
            field_source_table['zero_point']  = -999.0
            field_source_table['zpoint_err']  = -999.0
            mcd.output_log_entry(path_logfile,'Converting pixel coordinates done.')
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: get_radec_from_pixel_coords_array()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - RA/Dec computation from pixel coordinate array using get_radec_from_pixel_coords_array() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - RA/Dec computation from pixel coordinate array using get_radec_from_pixel_coords_array() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'RA/Dec computation from pixel coordinate array failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: get_radec_from_pixel_coords_array()')
    return field_source_table,proc_status,keep_going


##### FUNCTION DEFINITIONS -- SOURCE EXTRACTION #####

def measure_skylevel_stddev(fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    sky_mean,sky_stddev = 0.0,0.0
    if keep_going:
        try:
            image_data = fits.getdata(fits_filename_wcs)
            data_min = np.min(image_data)
            data_median = np.median(image_data)
            data_max = data_median + (data_median - data_min)
            data_array = image_data.flatten()
            data_array.sort(axis=0)
            num_pixels = len(data_array)
            idx,max_pix_value = 0,0
            while idx < num_pixels and max_pix_value < data_max:
                max_pix_value = data_array[idx]
                idx += 1
            num_pix_clipped = idx-1
            data_clipped = [0 for idx in range(num_pix_clipped)]
            for idx in range(0,num_pix_clipped):
                data_clipped[idx] = data_array[idx]
            sky_mean   = np.mean(data_clipped)
            sky_stddev = np.std(data_clipped)
            mcd.output_log_entry(path_logfile,'Sky level measured: mean={:.10f}, stddev={:.10f}'.format(sky_mean,sky_stddev))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: measure_skylevel_stddev()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Sky level measurement using measure_skylevel_stddev() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Sky level measurement using measure_skylevel_stddev() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Sky level measurement using measure_skylevel_stddev() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: measure_skylevel_stddev()')
    return sky_mean,sky_stddev,proc_status,keep_going


def extract_sources_sdss(sky_stddev,fits_filename_wcs,fits_filestem,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Run tphot to extract sources using both 2D Waussian fits and trailed Waussian fits
    output_path_stars,output_path_trails,output_path_moments,output_path_source_list = '','','',''
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Running tphot to extract sources from {:s} ...'.format(fits_filename_wcs))
            output_path_stars            = fits_filestem + '_wcs.stars'
            output_path_trails           = fits_filestem + '_wcs.trails'
            output_path_moments          = fits_filestem + '_wcs.moments'
            output_path_source_list      = fits_filestem + '_wcs.srclist'
            tphot_cmd = ''
            if os.path.isfile('/Users/hhsieh/Astro/tools/tphot/tphot'):
                tphot_cmd             = '/Users/hhsieh/Astro/tools/tphot/tphot'
            elif os.path.isfile('/home/hhsieh/tools/tphot/tphot'):
                tphot_cmd             = '/home/hhsieh/tools/tphot/tphot'  # on pohaku
            elif os.path.isfile('/atlas/bin/tphot'):
                tphot_cmd             = '/atlas/bin/tphot'  # on atlas
            cmd_bias              = '-1000'  # needed for SDSS because skymean ~ 0
            cmd_border            = '25'
            cmd_net               = '{:.3f}'.format(sky_stddev)
            cmd_min               = '{:.3f}'.format(sky_stddev)
            cmd_snr               = '3'
            cmd_move              = '10'
            cmd_okfit             = '0'

            # extract sources using two-dimensional Waussian fits
            cmd = [tphot_cmd,fits_filename_wcs,'-snr',cmd_snr,'-bias',cmd_bias,'-border',cmd_border,'-net',cmd_net,'-min',cmd_min,'-move',cmd_move,'-okfit',cmd_okfit,'-out','temp.out']
            mcd.output_log_entry(path_logfile,'{:s} {:s} -snr {:s} -bias {:s} -border {:s} -net {:s} -min {:s} -move {:s} -okfit {:s} -out temp.out'.format(tphot_cmd,fits_filename_wcs,cmd_snr,cmd_bias,cmd_border,cmd_net,cmd_min,cmd_move,cmd_okfit))
            process = subprocess.call(cmd)
            
            cmd = ['sort','-g','-r','-k','3','temp.out']
            mcd.output_log_entry(path_logfile,'sort -g -r -k 3 temp.out')
            with open('temp.srt','w') as of:
                process = subprocess.call(cmd,stdout=of)
            
            cmd = [tphot_cmd,fits_filename_wcs,'-obj','temp.srt','-out',output_path_stars,'-moment',output_path_moments,'-subtract','-snr',cmd_snr,'-bias',cmd_bias,'-border',cmd_border,'-net',cmd_net,'-min',cmd_min,'-move',cmd_move,'-okfit',cmd_okfit]
            mcd.output_log_entry(path_logfile,'{:s} {:s} -obj temp.srt -out {:s} -moment {:s} -subtract -snr {:s} -bias {:s} -border {:s} -net {:s} -min {:s} -move {:s} -okfit {:s}'.format(tphot_cmd,fits_filename_wcs,output_path_stars,output_path_moments,cmd_snr,cmd_bias,cmd_border,cmd_net,cmd_min,cmd_move,cmd_okfit))
            process = subprocess.call(cmd)
            
            os.remove('temp.out')
            os.remove('temp.srt')
            
            source_count = len(open(output_path_stars).readlines())  # count number of lines in file
            if source_count > 0:
                mcd.output_log_entry(path_logfile,'{:d} extracted sources using 2D Waussian fits written to {:s}'.format(source_count,output_path_stars))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'No sources extracted from {:s} using 2D Waussian fits.'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Source extraction using 2D Waussian fits using extract_sources_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Source extraction using 2D Waussian fits using extract_sources_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'No sources extracted using 2D Waussian fits'
                keep_going = False
            
            # extract source list
            x = np.genfromtxt(output_path_stars,skip_header=1,usecols=(0))
            y = np.genfromtxt(output_path_stars,skip_header=1,usecols=(1))
            with open(output_path_source_list,'w') as of:
                num_sources = len(x)
                for idx in range(0,num_sources):
                    of.write('{:7.2f}  {:7.2f}\n'.format(x[idx],y[idx]))

            # extract sources using trailed Waussian fits
            cmd = [tphot_cmd,fits_filename_wcs,'-trail','-obj',output_path_source_list,'-out',output_path_trails,'-subtract','-snr',cmd_snr,'-bias',cmd_bias,'-border',cmd_border,'-net',cmd_net,'-min',cmd_min,'-move',cmd_move,'-okfit',cmd_okfit,'-chin','10000']
            mcd.output_log_entry(path_logfile,'{:s} {:s} -trail -obj {:s} -out {:s} -subtract -snr {:s} -bias {:s} -border {:s} -net {:s} -min {:s} -move {:s} -okfit {:s} -chin 10000'.format(tphot_cmd,fits_filename_wcs,output_path_source_list,output_path_trails,cmd_snr,cmd_bias,cmd_border,cmd_net,cmd_min,cmd_move,cmd_okfit))
            process = subprocess.call(cmd)
            source_count = len(open(output_path_trails).readlines()) - 1  # count number of lines in file
            if source_count > 0:
                mcd.output_log_entry(path_logfile,'{:d} extracted sources using trailed Waussian fits written to {:s}'.format(source_count,output_path_trails))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'No sources extracted from {:s} using trailed Waussian fits.'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Source extraction using trailed Waussian fits using extract_sources_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Source extraction using trailed Waussian fits using extract_sources_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'No sources extracted using trailed Waussian fits'
                keep_going = False
                
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: extract_sources_sdss()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Source extraction using extract_sources_sdss() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Source extraction using extract_sources_sdss() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Source extraction using extract_sources_sdss() failed'
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: extract_sources_sdss()')
    return output_path_stars,output_path_trails,output_path_moments,output_path_source_list,proc_status,keep_going
    
    
def create_field_source_table(output_path_stars,output_path_trails,output_path_moments,fits_filename_wcs,match_threshold_pix,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Create calibration star tables
    field_source_table = Table()
    calib_source_table = Table()
    if keep_going:
        try:
            # Create field source table with trailed Waussian fits
            mcd.output_log_entry(path_logfile,'Creating field source table from {:s}...'.format(output_path_trails))
            x_t     = np.genfromtxt(output_path_trails,skip_header=1,usecols=(0))
            y_t     = np.genfromtxt(output_path_trails,skip_header=1,usecols=(1))
            flux_t  = np.genfromtxt(output_path_trails,skip_header=1,usecols=(7))
            dflux_t = np.genfromtxt(output_path_trails,skip_header=1,usecols=(8))
            trail_t = np.genfromtxt(output_path_trails,skip_header=1,usecols=(9))
            fwhm_t  = np.genfromtxt(output_path_trails,skip_header=1,usecols=(10))
            phi_t   = np.genfromtxt(output_path_trails,skip_header=1,usecols=(11))
            chiN_t  = np.genfromtxt(output_path_trails,skip_header=1,usecols=(13))
            field_source_table = Table([x_t,y_t,flux_t,dflux_t,trail_t,fwhm_t,phi_t,chiN_t], \
                                       names=('x_t','y_t','flux_t','dflux_t','trail_t','fwhm_t','phi_t','chiN_t'), \
                                       dtype=('f8','f8','f8','f8','f8','f8','f8','f8'))
            field_source_table['SNR_t']   = -999.0
            field_source_table['x_s']     = -999.0
            field_source_table['y_s']     = -999.0
            field_source_table['flux_s']  = -999.0
            field_source_table['dflux_s'] = -999.0
            field_source_table['major_s'] = -999.0
            field_source_table['minor_s'] = -999.0
            field_source_table['fwhm_s']  = -999.0
            field_source_table['phi_s']   = -999.0
            field_source_table['chiN_s']  = -999.0
            field_source_table['SNR_s']   = -999.0
            field_source_table['rKron_s'] = -999.0
            num_sources_trails = len(field_source_table)
            for idx in range(0,num_sources_trails):
                field_source_table[idx]['SNR_t'] = field_source_table[idx]['flux_t'] / field_source_table[idx]['dflux_t']
            mcd.output_log_entry(path_logfile,'Creating field source table done.')

            # Read in data for 2D Waussian fit results
            mcd.output_log_entry(path_logfile,'Reading in data from {:s} and {:s}...'.format(output_path_stars,output_path_moments))
            x_s     = np.genfromtxt(output_path_stars,skip_header=1,usecols=(0))
            y_s     = np.genfromtxt(output_path_stars,skip_header=1,usecols=(1))
            flux_s  = np.genfromtxt(output_path_stars,skip_header=1,usecols=(7))
            dflux_s = np.genfromtxt(output_path_stars,skip_header=1,usecols=(8))
            major_s = np.genfromtxt(output_path_stars,skip_header=1,usecols=(9))
            minor_s = np.genfromtxt(output_path_stars,skip_header=1,usecols=(10))
            phi_s   = np.genfromtxt(output_path_stars,skip_header=1,usecols=(11))
            chiN_s  = np.genfromtxt(output_path_stars,skip_header=1,usecols=(13))
            rKron_s = np.genfromtxt(output_path_moments,skip_header=1,usecols=(6))
            num_sources_stars = len(x_s)
            SNR_s   = [0 for idx in range(num_sources_stars)]
            for idx in range(0,num_sources_stars):
                SNR_s[idx] = flux_s[idx] / dflux_s[idx]
            mcd.output_log_entry(path_logfile,'Reading in data done.')
            
            # Add 2D Waussian fit data to field source table
            mcd.output_log_entry(path_logfile,'Adding data from 2D Waussian fit files to field source table...')
            for idx_t in range(num_sources_trails):
                idx_s = 0
                x_star,y_star,flux_star,dflux_star,major_star,minor_star,fwhm_star,phi_star,chiN_star,rKron_star,SNR_star = 0,0,0,0,0,0,0,0,0,0,0
                star_source_found = False
                while not star_source_found and idx_s < num_sources_stars:
                    if abs(x_s[idx_s]-field_source_table[idx_t]['x_t']) < match_threshold_pix and abs(y_s[idx_s]-field_source_table[idx_t]['y_t']) < match_threshold_pix:
                        x_star     = x_s[idx_s]
                        y_star     = y_s[idx_s]
                        flux_star  = flux_s[idx_s]
                        dflux_star = dflux_s[idx_s]
                        major_star = major_s[idx_s]
                        minor_star = minor_s[idx_s]
                        fwhm_star  = ((major_star**2 + minor_star**2)/2)**0.5
                        phi_star   = phi_s[idx_s]
                        chiN_star  = chiN_s[idx_s]
                        rKron_star = rKron_s[idx_s]
                        SNR_star   = SNR_s[idx_s]
                        star_source_found = True
                    else:
                        idx_s += 1
                field_source_table[idx_t]['x_s']     = x_star
                field_source_table[idx_t]['y_s']     = y_star
                field_source_table[idx_t]['flux_s']  = flux_star
                field_source_table[idx_t]['dflux_s'] = dflux_star
                field_source_table[idx_t]['major_s'] = major_star
                field_source_table[idx_t]['minor_s'] = minor_star
                field_source_table[idx_t]['fwhm_s']  = fwhm_star
                field_source_table[idx_t]['phi_s']   = phi_star
                field_source_table[idx_t]['chiN_s']  = chiN_star
                field_source_table[idx_t]['rKron_s'] = rKron_star
                field_source_table[idx_t]['SNR_s']   = SNR_star
            mcd.output_log_entry(path_logfile,'Adding data to field source table done.')
            
            # Create calibration source table
            calib_source_table = Table(names=('x_t','y_t','flux_t','dflux_t','trail_t','fwhm_t','phi_t','chiN_t','SNR_t', \
                'x_s','y_s','flux_s','dflux_s','major_s','minor_s','fwhm_s','phi_s','chiN_s','SNR_s','rKron_s'))
            for idx in range(len(field_source_table)):
                calib_source_table.add_row((field_source_table[idx]['x_t'],field_source_table[idx]['y_t'], \
                    field_source_table[idx]['flux_t'],field_source_table[idx]['dflux_t'], \
                    field_source_table[idx]['trail_t'],field_source_table[idx]['fwhm_t'], \
                    field_source_table[idx]['phi_t'],field_source_table[idx]['chiN_t'],field_source_table[idx]['SNR_t'], \
                    field_source_table[idx]['x_s'],field_source_table[idx]['y_s'], \
                    field_source_table[idx]['flux_s'],field_source_table[idx]['dflux_s'], \
                    field_source_table[idx]['major_s'],field_source_table[idx]['minor_s'],field_source_table[idx]['fwhm_s'], \
                    field_source_table[idx]['phi_s'],field_source_table[idx]['chiN_s'],field_source_table[idx]['SNR_s'], \
                    field_source_table[idx]['rKron_s']))

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: create_field_source_table'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Field source table creation using create_field_source_table() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Field source table creation using create_field_source_table() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Field source table creation using create_field_source_table() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_field_source_table')
    return field_source_table,calib_source_table,proc_status,keep_going


def remove_bad_sources(field_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Remove bad sources from tphot source list (e.g., flux < 0 or SNR < 3)
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Removing flux<0 and SNR<3 sources from tphot source list...')
            idx,sources_removed = 0,0
            num_field_sources = len(field_source_table)
            while idx < len(field_source_table):
                if field_source_table[idx]['flux_s'] < 0 or field_source_table[idx]['flux_t'] < 0 or field_source_table[idx]['SNR_s'] < 3:
                    field_source_table.remove_row(idx)
                    sources_removed += 1
                else:
                    idx += 1
            mcd.output_log_entry(path_logfile,'{:d} negative-flux or SNR < 3 sources removed ({:d} remaining).'.format(sources_removed,(num_field_sources - sources_removed)))
            if (num_field_sources - sources_removed) < 3:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Processing stopped for {:s}: remove_bad_sources() - < 3 field sources remaining'.format(fits_filename_wcs))
                print(e)
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Bad source removal using remove_bad_sources() stopped - < 3 field sources remaining.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Bad source removal using remove_bad_sources() stopped - < 3 field sources remaining.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'Fewer than 3 field sources remaining after bad source removal'
                keep_going = False
            else:
                keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: remove_bad_sources()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Bad source removal using remove_bad_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Bad source removal using remove_bad_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Bad source removal using remove_bad_sources() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: remove_bad_sources()')
    return field_source_table,proc_status,keep_going


def remove_lowsnr_sources(calib_star_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Remove low-SNR detections from calibration star table
    if keep_going:
        mcd.output_log_entry(path_logfile,'Removing low-SNR detections from calibration star table...')
        try:
            idx,sources_removed = 0,0
            num_calib_stars = len(calib_star_table)
            while idx < len(calib_star_table):
                if calib_star_table[idx]['SNR_s'] < 5:
                    calib_star_table.remove_row(idx)
                    sources_removed += 1
                else:
                    idx += 1
            mcd.output_log_entry(path_logfile,'Removing low-SNR detections done.')
            mcd.output_log_entry(path_logfile,'{:d} calibration sources with SNR < 5 removed ({:d} remaining).'.format(sources_removed,(num_calib_stars - sources_removed)))
            if (num_calib_stars - sources_removed) < 3:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Processing stopped for {:s}: remove_lowsnr_sources() - < 3 calibration sources remaining'.format(fits_filename_wcs))
                print(e)
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Low SNR source removal using remove_lowsnr_sources() stopped - < 3 calibration sources remaining.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Low SNR source removal using remove_lowsnr_sources() stopped - < 3 calibration sources remaining.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'Fewer than 3 calibration sources remaining after low SNR source removal'
                keep_going = False
            else:
                keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: remove_lowsnr_sources()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Low SNR source removal using remove_lowsnr_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Low SNR source removal using remove_lowsnr_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Low SNR source removal using remove_lowsnr_sources() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: remove_lowsnr_sources()')
    return calib_star_table,proc_status,keep_going


def compute_exp_source_density(exp_npix_x,exp_npix_y,wcs_pixscale_x,wcs_pixscale_y,field_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Compute average source density per square arcminute for entire exposure
    source_density = 0
    if keep_going:
        try:
            num_sources         = len(field_source_table)
            image_size_sqarcmin = exp_npix_x*wcs_pixscale_x/60 * exp_npix_y*wcs_pixscale_y/60
            source_density      = num_sources / image_size_sqarcmin
            mcd.output_log_entry(path_logfile,'Exposure source density: {:.3f}'.format(source_density))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: compute_exp_source_density()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Source density computation using compute_exp_source_density() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Source density computation using compute_exp_source_density() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Source density computation using compute_exp_source_density() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_exp_source_density()')
    return source_density,proc_status,keep_going


##### FUNCTION DEFINITIONS -- PHOTOMETRIC CALIBRATION #####

def retrieve_refcat_sources(ra_deg,dec_deg,radius_arcmin,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve refcat sources
    refcat_table = 0
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Retrieving refcat sources within {:.1f} arcmin of RA={:.6f},Dec={:.6f}'.format(radius_arcmin,ra_deg,dec_deg))
            filename = fits_filename_wcs[:-9]+'_refcat_sources.dat'
            outf = open(filename,'w')
            with open(filename,'w') as outf:
                outf.write(' objRA       objDec      gmag   gerr  rmag   rerr  imag   ierr  zmag   zerr\n')
                radius_deg = radius_arcmin / 60
                cmd_ra  = '{:f}'.format(ra_deg)
                cmd_dec = '{:f}'.format(dec_deg)
                cmd_radius = '{:f}'.format(radius_deg)
                cmd_refcat = ''
                if os.path.exists('/users/hhsieh/astro/tools/refcat'):
                    cmd_refcat = '/users/hhsieh/astro/tools/refcat'
                    directories = '/users/hhsieh/astro/tools/refcat_catalog/00_m_16_binary,/users/hhsieh/astro/tools/refcat_catalog/16_m_17_binary,/users/hhsieh/astro/tools/refcat_catalog/17_m_18_binary,/users/hhsieh/astro/tools/refcat_catalog/18_m_19_binary,/users/hhsieh/astro/tools/refcat_catalog/19_m_20_binary'
                elif os.path.exists('/home/hhsieh/tools/refcat'):
                    cmd_refcat = '/home/hhsieh/tools/refcat'  # on pohaku
                    directories = '/data2/refcat_catalog/00_m_16_binary,/data2/refcat_catalog/16_m_17_binary,/data2/refcat_catalog/17_m_18_binary,/data2/refcat_catalog/18_m_19_binary,/data2/refcat_catalog/19_m_20_binary'
                elif os.path.exists('/atlas/bin/refcat'): 
                    cmd_refcat = '/atlas/bin/refcat'  # on atlas
                    directories = '/home/hhsieh/tools/refcat_catalog/00_m_16_binary,/home/hhsieh/tools/refcat_catalog/16_m_17_binary,/home/hhsieh/tools/refcat_catalog/17_m_18_binary,/home/hhsieh/tools/refcat_catalog/18_m_19_binary,/home/hhsieh/tools/refcat_catalog/19_m_20_binary'
                var_fields = 'RA,Dec,g,dg,r,dr,i,di,z,dz'
                cmd_dir,cmd_rad,cmd_var,cmd_nohdr,cmd_bin = '-dir','-rad','-var','-nohdr','-bin'
                cmd = [cmd_refcat,cmd_ra,cmd_dec,cmd_dir,directories,cmd_bin,cmd_rad,cmd_radius,cmd_var,var_fields,cmd_nohdr]
                subprocess.call(cmd,stdout=outf)
            mcd.output_log_entry(path_logfile,'refcat catalog sources written to {:s}.'.format(filename))
            
            # parse local file into astropy.table object
            objRA  = np.genfromtxt(filename,skip_header=1,usecols=(0))
            objDec = np.genfromtxt(filename,skip_header=1,usecols=(1))
            gmag   = np.genfromtxt(filename,skip_header=1,usecols=(2))
            gerr   = np.genfromtxt(filename,skip_header=1,usecols=(3))
            rmag   = np.genfromtxt(filename,skip_header=1,usecols=(4))
            rerr   = np.genfromtxt(filename,skip_header=1,usecols=(5))
            imag   = np.genfromtxt(filename,skip_header=1,usecols=(6))
            ierr   = np.genfromtxt(filename,skip_header=1,usecols=(7))
            zmag   = np.genfromtxt(filename,skip_header=1,usecols=(8))
            zerr   = np.genfromtxt(filename,skip_header=1,usecols=(9))
            refcat_table = Table([objRA,objDec,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr], \
                names=('right_ascension','declination','gp1_mag','gp1_err','rp1_mag','rp1_err','ip1_mag','ip1_err','zp1_mag','zp1_err'), \
                dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
            num_refcat_sources = len(refcat_table)
            if num_refcat_sources > 0:
                idx = 0
                # Remove sources that are likely to be saturated or have dummy values (i.e., -999)
                while idx < len(refcat_table):
                    if refcat_table[idx]['gp1_mag'] < 10 or refcat_table[idx]['rp1_mag'] < 10 or refcat_table[idx]['ip1_mag'] < 10 or refcat_table[idx]['zp1_mag'] < 10:
                        refcat_table.remove_row(idx)
                    else:
                        idx += 1
                num_refcat_sources = len(refcat_table)
            if num_refcat_sources > 0:
                mcd.output_log_entry(path_logfile,'{:d} refcat catalog sources retrieved.'.format(num_refcat_sources))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'No refcat calibration sources found for {:s}'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - No refcat calibration sources found using retrieve_refcat_sources().\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - No refcat calibration sources found using retrieve_refcat_sources().\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'No refcat calibration sources found using retrieve_refcat_sources()'
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: retrieve_refcat_sources()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Refcat source retrieval using retrieve_refcat_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Refcat source retrieval using retrieve_refcat_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Refcat source retrieval using retrieve_refcat_sources() failed'
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_refcat_sources()')
    return refcat_table,proc_status,keep_going


def match_refcat_sources(field_source_table,refcat_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Find matches for sources extracted by tphot in PS1 source catalog
    calib_star_table = 0
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Finding matches for sources extracted by tphot in refcat catalog...')
            calib_star_table = Table(names=('x_s','y_s','flux_s','dflux_s','SNR_s', \
                'RA','Dec','refcat_RA','refcat_Dec','gp1_mag','gp1_err','rp1_mag','rp1_err','ip1_mag', \
                'ip1_err','zp1_mag','zp1_err','g_sdss_mag','g_sdss_err','r_sdss_mag','r_sdss_err','i_sdss_mag', \
                'i_sdss_err','z_sdss_mag','z_sdss_err','B_mag','B_err','V_mag','V_err','Rc_mag','Rc_err','Ic_mag', \
                'Ic_err','target_dist','zero_point','zpoint_err'))
            num_field_sources,num_refcat_stars = len(field_source_table),len(refcat_table)
            stars_matched = 0
            for idx1 in range(0,num_field_sources):
                idx2,star_matched = 0,0
                while idx2 < num_refcat_stars and star_matched == 0:
                    distance = ((field_source_table[idx1]['RA'] - refcat_table[idx2]['right_ascension'])**2 +
                        (field_source_table[idx1]['Dec'] - refcat_table[idx2]['declination'])**2)**0.5
                    if distance < 0.000417:  # If target matches within 1.5 arcsec
                        star_matched = 1
                        copy_fieldsource_row_to_calibstar_table(field_source_table,calib_star_table,idx1,fits_filename_wcs,thread_idx,keep_going,path_logfile,path_errorfile)
                        calib_star_table[stars_matched]['refcat_RA']  = float(refcat_table[idx2]['right_ascension'])
                        calib_star_table[stars_matched]['refcat_Dec'] = float(refcat_table[idx2]['declination'])
                        calib_star_table[stars_matched]['gp1_mag'] = refcat_table[idx2]['gp1_mag']
                        calib_star_table[stars_matched]['gp1_err'] = refcat_table[idx2]['gp1_err']
                        calib_star_table[stars_matched]['rp1_mag'] = refcat_table[idx2]['rp1_mag']
                        calib_star_table[stars_matched]['rp1_err'] = refcat_table[idx2]['rp1_err']
                        calib_star_table[stars_matched]['ip1_mag'] = refcat_table[idx2]['ip1_mag']
                        calib_star_table[stars_matched]['ip1_err'] = refcat_table[idx2]['ip1_err']
                        calib_star_table[stars_matched]['zp1_mag'] = refcat_table[idx2]['zp1_mag']
                        calib_star_table[stars_matched]['zp1_err'] = refcat_table[idx2]['zp1_err']
                        stars_matched += 1
                    else:
                        idx2 += 1
            mcd.output_log_entry(path_logfile,'{:d} extracted sources matched to refcat catalog sources for {:s}.'.format(stars_matched,fits_filename_wcs))
            if stars_matched == 0:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'No refcat sources matched to image sources using match_refcat_sources() for {:s}.'.format(fits_filename_wcs))
                with open(fits_filename_wcs[:-9]+'.refcat_source_matching_failed','w') as f:
                    f.write('{:s} - No refcat sources matched to image sources using match_refcat_sources().\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'No refcat sources matched to image sources using match_refcat_sources()'
                keep_going = False
            elif stars_matched < 3:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Insufficient refcat sources (n<3) refcat calibration sources matched for {:s}.'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Insufficient refcat sources (n<3) matched to image sources using match_refcat_sources().\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Insufficient refcat sources (n<3) matched to image sources using match_refcat_sources().\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'Fewer than 3 refcat sources matched to image sources using match_refcat_sources()'
                keep_going = False
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: match_refcat_sources()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Refcat source matching using match_refcat_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Refcat source matching using match_refcat_sources() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Refcat source matching using match_refcat_sources() failed'
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: match_refcat_sources()')
    return calib_star_table,proc_status,keep_going


def copy_fieldsource_row_to_calibstar_table(field_source_table,calib_star_table,target_idx,fits_filename_wcs,thread_idx,keep_going,path_logfile,path_errorfile):
    # Copy source from refcat source list to calibration star list
    if keep_going:
        try:
            calib_star_table.add_row((field_source_table[target_idx]['x_s'],field_source_table[target_idx]['y_s'], \
                field_source_table[target_idx]['flux_s'],field_source_table[target_idx]['dflux_s'], \
                field_source_table[target_idx]['SNR_s'], \
                field_source_table[target_idx]['RA'],field_source_table[target_idx]['Dec'], \
                field_source_table[target_idx]['ps1_RA'],field_source_table[target_idx]['ps1_Dec'], \
                field_source_table[target_idx]['gp1_mag'],field_source_table[target_idx]['gp1_err'], \
                field_source_table[target_idx]['rp1_mag'],field_source_table[target_idx]['rp1_err'], \
                field_source_table[target_idx]['ip1_mag'],field_source_table[target_idx]['ip1_err'], \
                field_source_table[target_idx]['zp1_mag'],field_source_table[target_idx]['zp1_err'], \
                field_source_table[target_idx]['g_sdss_mag'],field_source_table[target_idx]['g_sdss_err'], \
                field_source_table[target_idx]['r_sdss_mag'],field_source_table[target_idx]['r_sdss_err'], \
                field_source_table[target_idx]['i_sdss_mag'],field_source_table[target_idx]['i_sdss_err'], \
                field_source_table[target_idx]['z_sdss_mag'],field_source_table[target_idx]['z_sdss_err'], \
                field_source_table[target_idx]['B_mag'],field_source_table[target_idx]['B_err'], \
                field_source_table[target_idx]['V_mag'],field_source_table[target_idx]['V_err'], \
                field_source_table[target_idx]['Rc_mag'],field_source_table[target_idx]['Rc_err'], \
                field_source_table[target_idx]['Ic_mag'],field_source_table[target_idx]['Ic_err'], \
                field_source_table[target_idx]['target_dist'], \
                field_source_table[target_idx]['zero_point'],field_source_table[target_idx]['zpoint_err']))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: copy_fieldsource_row_to_calibstar_table()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Calibration star table insertion using copy_fieldsource_row_to_calibstar_table() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Calibration star table insertion using copy_fieldsource_row_to_calibstar_table() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: copy_fieldsource_row_to_calibstar_table()')
    return calib_star_table,keep_going


def convert_ps1_mags(gp1_mag,gp1_err,rp1_mag,rp1_err,ip1_mag,ip1_err,zp1_mag,zp1_err,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Convert PS1 magnitudes to SDSS and Johnson-Cousins systems
    magnitudes = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    if keep_going:
        try:
            # Transform magnitudes
            g_ps2sdss_a0,g_ps2sdss_a1,g_ps2sdss_a2,g_ps2sdss_b0,g_ps2sdss_b1 =  0.013, 0.145, 0.019, 0.014, 0.162
            r_ps2sdss_a0,r_ps2sdss_a1,r_ps2sdss_a2,r_ps2sdss_b0,r_ps2sdss_b1 = -0.001, 0.004, 0.007,-0.001, 0.011
            i_ps2sdss_a0,i_ps2sdss_a1,i_ps2sdss_a2,i_ps2sdss_b0,i_ps2sdss_b1 = -0.005, 0.011, 0.010,-0.004, 0.020
            z_ps2sdss_a0,z_ps2sdss_a1,z_ps2sdss_a2,z_ps2sdss_b0,z_ps2sdss_b1 =  0.013,-0.039,-0.012, 0.013,-0.050
            B_ps2jc_a0,  B_ps2jc_a1,  B_ps2jc_a2,  B_ps2jc_b0,  B_ps2jc_b1   =  0.212, 0.556, 0.034, 0.213, 0.587
            V_ps2jc_a0,  V_ps2jc_a1,  V_ps2jc_a2,  V_ps2jc_b0,  V_ps2jc_b1   =  0.005, 0.462, 0.013, 0.006, 0.474
            Rc_ps2jc_a0, Rc_ps2jc_a1, Rc_ps2jc_a2, Rc_ps2jc_b0, Rc_ps2jc_b1  = -0.137,-0.108,-0.029,-0.138,-0.131
            Ic_ps2jc_a0, Ic_ps2jc_a1, Ic_ps2jc_a2, Ic_ps2jc_b0, Ic_ps2jc_b1  = -0.366,-0.136,-0.018,-0.367,-0.149
            # y = A0 + A1x + A2x^2 = B0+B1x - Tonry et al. (2012, ApJ, 750, 99)
            x = gp1_mag - rp1_mag
            g_sdss_mag1 = gp1_mag + g_ps2sdss_a0 + g_ps2sdss_a1*x + g_ps2sdss_a2*(x**2)
            r_sdss_mag1 = rp1_mag + r_ps2sdss_a0 + r_ps2sdss_a1*x + r_ps2sdss_a2*(x**2)
            i_sdss_mag1 = ip1_mag + i_ps2sdss_a0 + i_ps2sdss_a1*x + i_ps2sdss_a2*(x**2)
            z_sdss_mag1 = zp1_mag + z_ps2sdss_a0 + z_ps2sdss_a1*x + z_ps2sdss_a2*(x**2)
            g_sdss_mag2 = gp1_mag + g_ps2sdss_b0 + g_ps2sdss_b1*x
            r_sdss_mag2 = rp1_mag + r_ps2sdss_b0 + r_ps2sdss_b1*x
            i_sdss_mag2 = ip1_mag + i_ps2sdss_b0 + i_ps2sdss_b1*x
            z_sdss_mag2 = zp1_mag + z_ps2sdss_b0 + z_ps2sdss_b1*x
            B_mag1  = gp1_mag + B_ps2jc_a0  + B_ps2jc_a1*x  + B_ps2jc_a2*(x**2)
            V_mag1  = rp1_mag + V_ps2jc_a0  + V_ps2jc_a1*x  + V_ps2jc_a2*(x**2)
            Rc_mag1 = rp1_mag + Rc_ps2jc_a0 + Rc_ps2jc_a1*x + Rc_ps2jc_a2*(x**2)
            Ic_mag1 = ip1_mag + Ic_ps2jc_a0 + Ic_ps2jc_a1*x + Ic_ps2jc_a2*(x**2)
            B_mag2  = gp1_mag + B_ps2jc_b0  + B_ps2jc_b1*x
            V_mag2  = rp1_mag + V_ps2jc_b0  + V_ps2jc_b1*x
            Rc_mag2 = rp1_mag + Rc_ps2jc_b0 + Rc_ps2jc_b1*x
            Ic_mag2 = ip1_mag + Ic_ps2jc_b0 + Ic_ps2jc_b1*x
            g_sdss_mag,g_sdss_err = mcd.magavg([g_sdss_mag1,g_sdss_mag2],[gp1_err,gp1_err])
            r_sdss_mag,r_sdss_err = mcd.magavg([r_sdss_mag1,r_sdss_mag2],[rp1_err,rp1_err])
            i_sdss_mag,i_sdss_err = mcd.magavg([i_sdss_mag1,i_sdss_mag2],[ip1_err,ip1_err])
            z_sdss_mag,z_sdss_err = mcd.magavg([z_sdss_mag1,z_sdss_mag2],[zp1_err,zp1_err])
            B_mag, B_err  = mcd.magavg([B_mag1, B_mag2] ,[gp1_err,gp1_err])
            V_mag, V_err  = mcd.magavg([V_mag1, V_mag2] ,[rp1_err,rp1_err])
            Rc_mag,Rc_err = mcd.magavg([Rc_mag1,Rc_mag2],[rp1_err,rp1_err])
            Ic_mag,Ic_err = mcd.magavg([Ic_mag1,Ic_mag2],[ip1_err,ip1_err])
            magnitudes = [g_sdss_mag,g_sdss_err,r_sdss_mag,r_sdss_err,i_sdss_mag,i_sdss_err,z_sdss_mag,z_sdss_err, \
                B_mag,B_err,V_mag,V_err,Rc_mag,Rc_err,Ic_mag,Ic_err]
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: convert_ps1_mags()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - PS1 magnitude conversion using convert_ps1_mags() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - PS1 magnitude conversion using convert_ps1_mags() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'PS1 magnitude conversion using convert_ps1_mags() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: convert_ps1_mags()')
    return magnitudes,proc_status,keep_going
    

def compute_zero_point(fits_filename_wcs,filter_name,calib_star_table,exptime,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Compute zero points
    avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars = 0,0,0
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Computing zero points for {:s}...'.format(fits_filename_wcs))
            num_calib_stars = len(calib_star_table)
            for idx in range(0,num_calib_stars):
                magnitudes,proc_status,keep_going = convert_ps1_mags(calib_star_table[idx]['gp1_mag'],calib_star_table[idx]['gp1_err'], \
                    calib_star_table[idx]['rp1_mag'],calib_star_table[idx]['rp1_err'], \
                    calib_star_table[idx]['ip1_mag'],calib_star_table[idx]['ip1_err'], \
                    calib_star_table[idx]['zp1_mag'],calib_star_table[idx]['zp1_err'],fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                if keep_going:
                    keep_going = False
                    if filter_name == 'g':
                        calib_star_table[idx]['g_sdss_mag'] = magnitudes[0]
                        calib_star_table[idx]['g_sdss_err'] = magnitudes[1]
                        targ_magerr = ufloat(calib_star_table[idx]['g_sdss_mag'],calib_star_table[idx]['g_sdss_err'])
                    elif filter_name == 'r':
                        calib_star_table[idx]['r_sdss_mag'] = magnitudes[2]
                        calib_star_table[idx]['r_sdss_err'] = magnitudes[3]
                        targ_magerr = ufloat(calib_star_table[idx]['r_sdss_mag'],calib_star_table[idx]['r_sdss_err'])
                    elif filter_name == 'i':
                        calib_star_table[idx]['i_sdss_mag'] = magnitudes[4]
                        calib_star_table[idx]['i_sdss_err'] = magnitudes[5]
                        targ_magerr = ufloat(calib_star_table[idx]['i_sdss_mag'],calib_star_table[idx]['i_sdss_err'])
                    elif filter_name == 'z':
                        calib_star_table[idx]['z_sdss_mag'] = magnitudes[6]
                        calib_star_table[idx]['z_sdss_err'] = magnitudes[7]
                        targ_magerr = ufloat(calib_star_table[idx]['z_sdss_mag'],calib_star_table[idx]['z_sdss_err'])
                    elif filter_name == 'B':
                        calib_star_table[idx]['B_mag']      = magnitudes[8]
                        calib_star_table[idx]['B_err']      = magnitudes[9]
                        targ_magerr = ufloat(calib_star_table[idx]['B_mag'], calib_star_table[idx]['B_err'])
                    elif filter_name == 'V':
                        calib_star_table[idx]['V_mag']      = magnitudes[10]
                        calib_star_table[idx]['V_err']      = magnitudes[11]
                        targ_magerr = ufloat(calib_star_table[idx]['V_mag'], calib_star_table[idx]['V_err'])
                    elif filter_name == 'Rc':
                        calib_star_table[idx]['Rc_mag']     = magnitudes[12]
                        calib_star_table[idx]['Rc_err']     = magnitudes[13]
                        targ_magerr = ufloat(calib_star_table[idx]['Rc_mag'],calib_star_table[idx]['Rc_err'])
                    elif filter_name == 'Ic':
                        calib_star_table[idx]['Ic_mag']     = magnitudes[14]
                        calib_star_table[idx]['Ic_err']     = magnitudes[15]
                        targ_magerr = ufloat(calib_star_table[idx]['Ic_mag'],calib_star_table[idx]['Ic_err'])
                    star_fluxerr = ufloat(calib_star_table[idx]['flux_s'],calib_star_table[idx]['dflux_s'])
                    zero_point = targ_magerr + 2.5 * unumpy.log10(star_fluxerr/exptime)
                    calib_star_table[idx]['zero_point'] = zero_point.nominal_value
                    calib_star_table[idx]['zpoint_err'] = zero_point.std_dev
                    keep_going = True
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'Conversion of refcat source magnitudes to other filters failed for {:s}.'.format(fits_filename_wcs))
                    if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                        with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                            f.write('{:s} - Refcat source magnitude conversion using compute_zero_point() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                    else:
                        with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                            f.write('{:s} - Refcat source magnitude conversion using compute_zero_point() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                    proc_status = 'Refcat source magnitude conversion using compute_zero_point() failed'
                    keep_going = False

            if keep_going:
                mcd.output_log_entry(path_logfile,'Removing sources with outlier zero point magnitudes...')
                zeropoints = calib_star_table['zero_point']
                mu,sigma,proc_status,keep_going = generate_zeropoint_histogram(calib_star_table['zero_point'],fits_filename_wcs,'histogram1',proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                zeropoint_median = statistics.median(zeropoints)
                idx,zpoints_removed = 0,0
                num_calib_stars = len(calib_star_table)                
                while idx < len(calib_star_table):
                    # Remove stars with zero points more than 0.3 mag above or below the median
                    if calib_star_table[idx]['zero_point'] > (zeropoint_median+0.3) or calib_star_table[idx]['zero_point'] < (zeropoint_median-0.3):
                        calib_star_table.remove_row(idx)
                        zpoints_removed += 1
                    else:
                        idx += 1
                mcd.output_log_entry(path_logfile,'{:d} outlier zeropoints removed, {:d} zeropoints remaining.'.format(zpoints_removed,(num_calib_stars-zpoints_removed)))
                if (num_calib_stars-zpoints_removed) < 3:
                    mcd.output_log_error_entry(path_logfile,path_errorfile,'Insufficient zeropoints (n<3) remaining for compute_zero_point() to continue for {:s}'.format(fits_filename_wcs))
                    if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                        with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                            f.write('{:s} - Insufficient zeropoints (n<3) remaining for compute_zero_point() to continue.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                    else:
                        with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                            f.write('{:s} - Insufficient zeropoints (n<3) remaining for compute_zero_point() to continue.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                    proc_status = 'Insufficient zeropoints (n<3) remaining for compute_zero_point() to continue'
                    keep_going = False
                else:
                    mu,sigma,proc_status,keep_going = generate_zeropoint_histogram(calib_star_table['zero_point'],fits_filename_wcs,'histogram2',proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                    avg_zeropoint_mag,avg_zeropoint_err = mcd.magavg(calib_star_table['zero_point'],calib_star_table['zpoint_err'])
                    mcd.output_log_entry(path_logfile,'Number of stars used for photometric solution: {:d}'.format(num_calib_stars))
                    mcd.output_log_entry(path_logfile,'Zero point median: {:.3f}'.format(zeropoint_median))
                    mcd.output_log_entry(path_logfile,'Zero point magnitude solution: {:.3f}+/-{:.3f}'.format(avg_zeropoint_mag,avg_zeropoint_err))
                    with open(fits_filename_wcs[:-9]+'.photsolved','w') as f:
                        f.write('{:s} - Photometric calibration successful.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                    proc_status = 'Photometric calibration successful'
                    keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Zero point computation failed for {:s}.'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Zero point computation failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Zero point computation failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                proc_status = 'Zero point computation failed'
                keep_going = False
        except Exception as e:
            mcd.output_log_error_entry(path_logfile,path_errorfile,'Function failed for {:s}: compute_zero_point()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Zero point computation failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Zero point computation failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Zero point computation failed'
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_zero_point()')
    return calib_star_table,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,proc_status,keep_going


def generate_zeropoint_histogram(zeropoints,fits_filename_wcs,filename_suffix,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    mu,sigma = 0,0
    if keep_going:
        try:
            filename = fits_filename_wcs[:-9]+'_zeropoint_'+filename_suffix
            rect_plot_1 = [0.12,0.12,0.40,0.30]
            plt.figure(1,figsize=(12,12))
            axPlot1 = plt.axes(rect_plot_1)
            n_zeropoint_bins = 20
            zeropoint_min = np.amin(zeropoints)
            zeropoint_max = np.amax(zeropoints)
            zeropoint_bins = np.arange(zeropoint_min,zeropoint_max,(zeropoint_max-zeropoint_min)/n_zeropoint_bins)
            mu,sigma = norm.fit(zeropoints)
            x1,bins1,p1 = axPlot1.hist(zeropoints,bins=zeropoint_bins,histtype='bar',density=1,edgecolor='black',color='#2077b4',zorder=1)
            #y = mlab.normpdf(zeropoint_bins,mu,sigma)
            y = scipy.stats.norm.pdf(zeropoint_bins,mu,sigma)
            axPlot1.plot(zeropoint_bins,y,ls='--',lw=2,color='#880000',zorder=2)
            plt.savefig(filename + '.pdf',format='pdf',transparent=True)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Zero-point histogram generation and Gaussian fitting (mu={:.2f} mag; sigma={:.2f}) completed successfully.'.format(mu,sigma))
            keep_going = True
        except Exception as e:
            mcd.output_log_error_entry(path_logfile,path_errorfile,'Function failed for {:s}: generate_zeropoint_histogram()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Zero point histogram generation using generate_zeropoint_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Zero point histogram generation using generate_zeropoint_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Zero point histogram generation using generate_zeropoint_histogram() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: generate_zeropoint_histogram()')
    return mu,sigma,proc_status,keep_going


def compute_mags_from_zpoint(target_flux,target_fluxerr,exptime,zeropoint_mag,zeropoint_magerr,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Compute magnitude of source from zero point
    target_mag,target_err = 0,0
    if keep_going:
        try:
            target_flux_ufloat = ufloat(target_flux,target_fluxerr)
            zpoint_mag_ufloat  = ufloat(zeropoint_mag,zeropoint_magerr)
            target_mag_ufloat  = zpoint_mag_ufloat - 2.5 * unumpy.log10(target_flux_ufloat/exptime)
            target_mag         = target_mag_ufloat.nominal_value
            target_err         = target_mag_ufloat.std_dev
            keep_going = True
        except Exception as e:
            mcd.output_log_error_entry(path_logfile,path_errorfile,'Function failed: compute_mags_from_zpoint() (Flux={:f}+/-{:f}, ExpTime={:f}, ZPt={:f}+/-{:f})'.format(target_flux,target_fluxerr,exptime,zeropoint_mag,zeropoint_magerr))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Magnitude computation from zero point using compute_mags_from_zpoint() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Magnitude computation from zero point using compute_mags_from_zpoint() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Magnitude computation from zero point using compute_mags_from_zpoint() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_mags_from_zpoint()')
    return target_mag,target_err,proc_status,keep_going


def compute_detection_limits(sky_stddev,psf_width_mean_pix,wcs_pixscale_x,exptime,zeropoint_mag,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    limit_mag_ps,limit_mag_sb = 0,0
    if keep_going:
        try:
            ps_apert_noise = (math.pi*(2.5*psf_width_mean_pix)**2)**0.5 * sky_stddev  # using circular aperture with r=2.5*psf_width_mean_pix
            sb_apert_noise = ((1/wcs_pixscale_x)**2)**0.5 * sky_stddev                # using 1 arcsec x 1 arcsec aperture
            limit_flux_ps  = ps_apert_noise * 3
            limit_flux_sb  = sb_apert_noise * 3
            limit_mag_ps   = zeropoint_mag - 2.5*np.log10(limit_flux_ps/exptime)
            limit_mag_sb   = zeropoint_mag - 2.5*np.log10(limit_flux_sb/exptime)
            mcd.output_log_entry(path_logfile,'Point source detection limit:       {:.3f} mag'.format(limit_mag_ps))
            mcd.output_log_entry(path_logfile,'Surface brightness detection limit: {:.3f} mag/arcsec^2'.format(limit_mag_sb))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: compute_detection_limits()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Detection limits computation using compute_detection_limits() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Detection limits computation using compute_detection_limits() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Detection limits computation using compute_detection_limits() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_detection_limits()')
    return limit_mag_ps,limit_mag_sb,proc_status,keep_going


##### FUNCTION DEFINITIONS -- MULTI-APERTURE PHOTOMETRY #####

def extract_sources_multiap_tphot_sdss(fits_filename_wcs,apfactor,rKron_mean_pix,output_path,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile):
    tphot_cmd = ''
    if os.path.isfile('/Users/hhsieh/Astro/tools/tphot/tphot'):
        tphot_cmd             = '/Users/hhsieh/Astro/tools/tphot/tphot'
    elif os.path.isfile('/home/hhsieh/tools/tphot/tphot'):
        tphot_cmd             = '/home/hhsieh/tools/tphot/tphot'  # on pohaku
    elif os.path.isfile('/atlas/bin/tphot'):
        tphot_cmd             = '/atlas/bin/tphot'  # on atlas
    value_bias              = '-1000'  # needed for SDSS because skymean ~ 0
    value_border            = '25'
    value_net               = '{:.3f}'.format(sky_stddev)
    value_min               = '{:.3f}'.format(sky_stddev)
    value_snr               = '3'
    value_move              = '10'
    value_okfit             = '0'
    aprad = apfactor*rKron_mean_pix
    if (aprad + 20) < 40: skyrad = 40
    else: skyrad = aprad+20
    cmd = [tphot_cmd,fits_filename_wcs,'-trail','-obj',output_path_source_list,'-out',output_path,'-subtract','-snr',value_snr,'-aprad','{:.3f}'.format(aprad),'-bias',value_bias,'-border',value_border,'-net',value_net,'-min',value_min,'-move',value_move,'-okfit',value_okfit,'-chin','10000']
    mcd.output_log_entry(path_logfile,'{:s} {:s} -trail -obj {:s} -out {:s} -subtract -snr {:s} -aprad {:.3f} -bias {:s} -border {:s} -net {:s} -min {:s} -move {:s} -okfit {:s} -chin 10000'.format(tphot_cmd,fits_filename_wcs,output_path_source_list,output_path,value_snr,aprad,value_bias,value_border,value_net,value_min,value_move,value_okfit))
    process = subprocess.call(cmd)
    source_count = len(open(output_path).readlines()) - 1  # count number of lines in file
    if source_count > 0:
        mcd.output_log_entry(path_logfile,'{:d} extracted sources using trailed Waussian fits with aperture radius of {:.1f}*rKron written to {:s} .'.format(source_count,apfactor,output_path))
    else:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'No sources extracted from {:s} using trailed Waussian fits with aperture radius of {:.1f}*rKron.'.format(fits_filename_wcs,apfactor))
        if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
            with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                f.write('{:s} - Source extraction using trailed Waussian fits with aperture radius of {:.1f}*rKron using extract_sources_multiap_tphot_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),apfactor))
        else:
            with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                f.write('{:s} - Source extraction using trailed Waussian fits with aperture radius of {:.1f}*rKron using extract_sources_multiap_tphot_sdss() failed (no sources extracted).\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),apfactor))
    return None

def perform_multiap_photometry_tphot(instrument,field_source_table,rKron_mean_pix,sky_stddev,output_path_source_list,fits_filename_wcs,fits_filestem,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Run tphot to extract sources using both 2D Waussian fits and trailed Waussian fits using multiple apertures
    multiap_photom_table = ''
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Running tphot to extract sources from {:s}...'.format(fits_filename_wcs))
            output_path_05rKron   = fits_filestem + '_wcs.05rKron.trails'
            output_path_10rKron   = fits_filestem + '_wcs.10rKron.trails'
            output_path_20rKron   = fits_filestem + '_wcs.20rKron.trails'
            output_path_30rKron   = fits_filestem + '_wcs.30rKron.trails'
            output_path_40rKron   = fits_filestem + '_wcs.40rKron.trails'
            output_path_50rKron   = fits_filestem + '_wcs.50rKron.trails'
            tphot_cmd = ''
            if os.path.isfile('/Users/hhsieh/Astro/tools/tphot/tphot'):
                tphot_cmd             = '/Users/hhsieh/Astro/tools/tphot/tphot'
            elif os.path.isfile('/home/hhsieh/tools/tphot/tphot'):
                tphot_cmd             = '/home/hhsieh/tools/tphot/tphot'  # on pohaku
            elif os.path.isfile('/atlas/bin/tphot'):
                tphot_cmd             = '/atlas/bin/tphot'  # on atlas

            if instrument=='SDSS':
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,0.5,rKron_mean_pix,output_path_05rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,1.0,rKron_mean_pix,output_path_10rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,2.0,rKron_mean_pix,output_path_20rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,3.0,rKron_mean_pix,output_path_30rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,4.0,rKron_mean_pix,output_path_40rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
                extract_sources_multiap_tphot_sdss(fits_filename_wcs,5.0,rKron_mean_pix,output_path_50rKron,output_path_source_list,tphot_cmd,sky_stddev,path_logfile,path_errorfile)
            
            # read in multi-aperture photometry data from files
            x_05rK     = np.genfromtxt(output_path_05rKron,skip_header=1,usecols=(0))
            y_05rK     = np.genfromtxt(output_path_05rKron,skip_header=1,usecols=(1))
            flux_05rK  = np.genfromtxt(output_path_05rKron,skip_header=1,usecols=(7))
            dflux_05rK = np.genfromtxt(output_path_05rKron,skip_header=1,usecols=(8))
            x_10rK     = np.genfromtxt(output_path_10rKron,skip_header=1,usecols=(0))
            y_10rK     = np.genfromtxt(output_path_10rKron,skip_header=1,usecols=(1))
            flux_10rK  = np.genfromtxt(output_path_10rKron,skip_header=1,usecols=(7))
            dflux_10rK = np.genfromtxt(output_path_10rKron,skip_header=1,usecols=(8))
            x_20rK     = np.genfromtxt(output_path_20rKron,skip_header=1,usecols=(0))
            y_20rK     = np.genfromtxt(output_path_20rKron,skip_header=1,usecols=(1))
            flux_20rK  = np.genfromtxt(output_path_20rKron,skip_header=1,usecols=(7))
            dflux_20rK = np.genfromtxt(output_path_20rKron,skip_header=1,usecols=(8))
            x_30rK     = np.genfromtxt(output_path_30rKron,skip_header=1,usecols=(0))
            y_30rK     = np.genfromtxt(output_path_30rKron,skip_header=1,usecols=(1))
            flux_30rK  = np.genfromtxt(output_path_30rKron,skip_header=1,usecols=(7))
            dflux_30rK = np.genfromtxt(output_path_30rKron,skip_header=1,usecols=(8))
            x_40rK     = np.genfromtxt(output_path_40rKron,skip_header=1,usecols=(0))
            y_40rK     = np.genfromtxt(output_path_40rKron,skip_header=1,usecols=(1))
            flux_40rK  = np.genfromtxt(output_path_40rKron,skip_header=1,usecols=(7))
            dflux_40rK = np.genfromtxt(output_path_40rKron,skip_header=1,usecols=(8))
            x_50rK     = np.genfromtxt(output_path_50rKron,skip_header=1,usecols=(0))
            y_50rK     = np.genfromtxt(output_path_50rKron,skip_header=1,usecols=(1))
            flux_50rK  = np.genfromtxt(output_path_50rKron,skip_header=1,usecols=(7))
            dflux_50rK = np.genfromtxt(output_path_50rKron,skip_header=1,usecols=(8))

            # add multi-aperture photometry data to table
            multiap_photom_table = Table(names=('x_t','y_t', \
                'x_05rK','y_05rK','flux_t_05rK','dflux_t_05rK', \
                'x_10rK','y_10rK','flux_t_10rK','dflux_t_10rK', \
                'x_20rK','y_20rK','flux_t_20rK','dflux_t_20rK', \
                'x_30rK','y_30rK','flux_t_30rK','dflux_t_30rK', \
                'x_40rK','y_40rK','flux_t_40rK','dflux_t_40rK', \
                'x_50rK','y_50rK','flux_t_50rK','dflux_t_50rK'))
            num_field_sources = len(field_source_table)
            for idx_fs in range(0,num_field_sources):  # go through field source table and match to sources in multi-ap photometry files
                flux_trail_05rK,flux_trail_10rK,flux_trail_20rK,flux_trail_30rK,flux_trail_40rK,flux_trail_50rK       = 0,0,0,0,0,0
                dflux_trail_05rK,dflux_trail_10rK,dflux_trail_20rK,dflux_trail_30rK,dflux_trail_40rK,dflux_trail_50rK = 0,0,0,0,0,0
                x_fs   = field_source_table[idx_fs]['x_t']
                y_fs   = field_source_table[idx_fs]['y_t']
                
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_05rK) and not source_found:
                    if abs(x_fs - x_05rK[idx_ph]) < 1 and abs(y_fs - y_05rK[idx_ph]) < 1:
                        x_trail_05rK     = x_05rK[idx_ph]
                        y_trail_05rK     = y_05rK[idx_ph]
                        flux_trail_05rK  = flux_05rK[idx_ph]
                        dflux_trail_05rK = dflux_05rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1
                        
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_10rK) and not source_found:
                    if abs(x_fs - x_10rK[idx_ph]) < 1 and abs(y_fs - y_10rK[idx_ph]) < 1:
                        x_trail_10rK     = x_10rK[idx_ph]
                        y_trail_10rK     = y_10rK[idx_ph]
                        flux_trail_10rK  = flux_10rK[idx_ph]
                        dflux_trail_10rK = dflux_10rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1
                        
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_20rK) and not source_found:
                    if abs(x_fs - x_20rK[idx_ph]) < 1 and abs(y_fs - y_20rK[idx_ph]) < 1:
                        x_trail_20rK     = x_20rK[idx_ph]
                        y_trail_20rK     = y_20rK[idx_ph]
                        flux_trail_20rK  = flux_20rK[idx_ph]
                        dflux_trail_20rK = dflux_20rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1
                        
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_30rK) and not source_found:
                    if abs(x_fs - x_30rK[idx_ph]) < 1 and abs(y_fs - y_30rK[idx_ph]) < 1:
                        x_trail_30rK     = x_30rK[idx_ph]
                        y_trail_30rK     = y_30rK[idx_ph]
                        flux_trail_30rK  = flux_30rK[idx_ph]
                        dflux_trail_30rK = dflux_30rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1
                        
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_40rK) and not source_found:
                    if abs(x_fs - x_40rK[idx_ph]) < 1 and abs(y_fs - y_40rK[idx_ph]) < 1:
                        x_trail_40rK     = x_40rK[idx_ph]
                        y_trail_40rK     = y_40rK[idx_ph]
                        flux_trail_40rK  = flux_40rK[idx_ph]
                        dflux_trail_40rK = dflux_40rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1
                        
                idx_ph = 0
                source_found = False
                while idx_ph < len(flux_50rK) and not source_found:
                    if abs(x_fs - x_50rK[idx_ph]) < 1 and abs(y_fs - y_50rK[idx_ph]) < 1:
                        x_trail_50rK     = x_50rK[idx_ph]
                        y_trail_50rK     = y_50rK[idx_ph]
                        flux_trail_50rK  = flux_50rK[idx_ph]
                        dflux_trail_50rK = dflux_50rK[idx_ph]
                        source_found = True
                    else:
                        idx_ph += 1

                multiap_photom_table.add_row((x_fs,y_fs,x_trail_05rK,y_trail_05rK,flux_trail_05rK,dflux_trail_05rK, \
                    x_trail_10rK,y_trail_10rK,flux_trail_10rK,dflux_trail_10rK, \
                    x_trail_20rK,y_trail_20rK,flux_trail_20rK,dflux_trail_20rK, \
                    x_trail_30rK,y_trail_30rK,flux_trail_30rK,dflux_trail_30rK, \
                    x_trail_40rK,y_trail_40rK,flux_trail_40rK,dflux_trail_40rK, \
                    x_trail_50rK,y_trail_50rK,flux_trail_50rK,dflux_trail_50rK))
                
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: perform_multiap_photometry_tphot()'.format(fits_filename_wcs))
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Multi-aperture photometry using perform_multiap_photometry_tphot() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Multi-aperture photometry using perform_multiap_photometry_tphot() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Multi-aperture photometry using perform_multiap_photometry_tphot() failed'
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: perform_multiap_photometry_tphot()')
    return multiap_photom_table,proc_status,keep_going


##### FUNCTION DEFINITIONS -- PSF ANALYSIS #####

def compute_psf_width_mean(wcs_pixscale_x,calib_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Remove sources with potentially extended PSFs
    mu,mu_pix = 0,0
    if keep_going:
        mcd.output_log_entry(path_logfile,'Removing sources with potentially extended PSFs for {:s}...'.format(fits_filename_wcs))
        try:
            sources_removed = 0
            fwhm = calib_source_table['fwhm_s']
            mu,mu_pix,sigma,keep_going = generate_fwhm_histogram(fwhm,wcs_pixscale_x,fits_filename_wcs,'histogram1',thread_idx,keep_going,path_logfile,path_errorfile)
            fwhm_median = statistics.median(fwhm)
            idx = 0
            num_calib_sources = len(calib_source_table)
            while idx < len(calib_source_table):
                if calib_source_table[idx]['fwhm_s'] > (1.2 * fwhm_median):
                    calib_source_table.remove_row(idx)
                    sources_removed += 1
                else:
                    idx += 1
            mcd.output_log_entry(path_logfile,'{:d} sources with potentially extended PSFs (1.1 x median) removed ({:d} remaining).'.format(sources_removed,(num_calib_sources - sources_removed)))
            if (num_calib_sources - sources_removed) < 3:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Processing stopped for {:s}: compute_psf_width_mean() - < 3 calibration sources remaining after removal of sources with extended PSFs'.format(fits_filename_wcs))
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Removal of sources with extended PSFs using compute_psf_width_mean() stopped - < 3 calibration sources remaining after removal of sources with extended PSFs.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Removal of sources with extended PSFs using compute_psf_width_mean() stopped - < 3 calibration sources remaining after removal of sources with extended PSFs.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                keep_going = False
            else:
                fwhm = calib_source_table['fwhm_s']
                mu,mu_pix,sigma,keep_going = generate_fwhm_histogram(fwhm,wcs_pixscale_x,fits_filename_wcs,'histogram2',thread_idx,keep_going,path_logfile,path_errorfile)
                mcd.output_log_entry(path_logfile,'Average PSF width (with outlier rejection): {:.3f} arcsec'.format(mu))
                keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: compute_psf_width_mean()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Mean PSF width computation using compute_psf_width_mean() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Mean PSF width computation using compute_psf_width_mean() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Mean PSF width computation using compute_psf_width_mean() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_psf_width_mean()')
    return mu,mu_pix,calib_source_table,proc_status,keep_going


def generate_fwhm_histogram(fwhms,wcs_pixscale_x,fits_filename_wcs,file_suffix,thread_idx,keep_going,path_logfile,path_errorfile):
    mu,sigma = 0.0,0.0
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Performing FWHM histogram generation and Gaussian fitting for {:s}...'.format(fits_filename_wcs))
            filename = fits_filename_wcs[:-9]+'_fwhm_'+file_suffix
            rect_plot_1 = [0.12,0.12,0.40,0.30]
            plt.figure(1,figsize=(12,12))
            axPlot1 = plt.axes(rect_plot_1)
            fwhms = fwhms*wcs_pixscale_x
            n_fwhm_bins = 100
            fwhm_min = 0.5
            fwhm_max = 2.5
            fwhm_bins = np.arange(fwhm_min,fwhm_max,(fwhm_max-fwhm_min)/n_fwhm_bins)
            mu,sigma = norm.fit(fwhms)
            mu_pix = mu / wcs_pixscale_x
            x1,bins1,p1 = axPlot1.hist(fwhms,bins=fwhm_bins,histtype='bar',density=1,edgecolor='black',color='#2077b4',zorder=1)
            y = scipy.stats.norm.pdf(fwhm_bins,mu,sigma)
            axPlot1.plot(fwhm_bins,y,ls='--',lw=2,color='#880000',zorder=2)
            plt.savefig(filename + '.pdf',format='pdf',transparent=True)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'FWHM histogram generation and Gaussian fitting (mu={:.2f} arcsec; sigma={:.2f}) for {:s} completed successfully.'.format(mu,sigma,fits_filename_wcs))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: generate_fwhm_histogram()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - FWHM histogram generation using generate_fwhm_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - FWHM histogram generation using generate_fwhm_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: generate_fwhm_histogram()')
    return mu,mu_pix,sigma,keep_going


def compute_rKron_mean(wcs_pixscale_x,calib_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    # Compute mean Kron radius
    mu,mu_pix = 0,0
    if keep_going:
        mcd.output_log_entry(path_logfile,'Removing sources with potentially extended Kron radii for {:s}...'.format(fits_filename_wcs))
        try:
            sources_removed = 0
            rKron = calib_source_table['rKron_s']
            mu,mu_pix,sigma,keep_going = generate_rKron_histogram(rKron,wcs_pixscale_x,fits_filename_wcs,'histogram1',thread_idx,keep_going,path_logfile,path_errorfile)
            rKron_median = statistics.median(rKron)
            idx = 0
            num_calib_sources = len(calib_source_table)
            while idx < len(calib_source_table):
                if calib_source_table[idx]['rKron_s'] > (1.2 * rKron_median):
                    calib_source_table.remove_row(idx)
                    sources_removed += 1
                else:
                    idx += 1
            mcd.output_log_entry(path_logfile,'{:d} sources with potentially extended Kron radii (1.1 x median) removed ({:d} remaining).'.format(sources_removed,(num_calib_sources - sources_removed)))
            if (num_calib_sources - sources_removed) < 3:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Processing stopped for {:s}: compute_rKron_mean() - < 3 calibration sources remaining after removal of sources with extended Kron radii.'.format(fits_filename_wcs))
                print(e)
                if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                    with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                        f.write('{:s} - Removal of sources with extended Kron radii using compute_rKron_mean() stopped - < 3 calibration sources remaining after removal of sources with extended Kron radii.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                else:
                    with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                        f.write('{:s} - Removal of sources with extended Kron radii using compute_rKron_mean() stopped - < 3 calibration sources remaining after removal of sources with extended Kron radii.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
                keep_going = False
            else:
                rKron = calib_source_table['rKron_s']
                mu,mu_pix,sigma,keep_going = generate_rKron_histogram(rKron,wcs_pixscale_x,fits_filename_wcs,'histogram2',thread_idx,keep_going,path_logfile,path_errorfile)
                mcd.output_log_entry(path_logfile,'Mean Kron radius (with outlier rejection): {:.3f} arcsec'.format(mu))
                keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: compute_rKron_mean()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Mean Kron radius computation using compute_rKron_mean() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Mean Kron radius computation using compute_rKron_mean() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Mean Kron radius computation using compute_rKron_mean() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: compute_rKron_mean()')
    return mu,mu_pix,calib_source_table,proc_status,keep_going


def generate_rKron_histogram(rKron,wcs_pixscale_x,fits_filename_wcs,file_suffix,thread_idx,keep_going,path_logfile,path_errorfile):
    mu,sigma = 0.0,0.0
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Performing rKron histogram generation and Gaussian fitting for {:s}.'.format(fits_filename_wcs))
            filename = fits_filename_wcs[:-9]+'_rKron_'+file_suffix
            rect_plot_1 = [0.12,0.12,0.40,0.30]
            plt.figure(1,figsize=(12,12))
            axPlot1 = plt.axes(rect_plot_1)
            rKron = rKron*wcs_pixscale_x
            n_rKron_bins = 100
            rKron_min = 0.5
            rKron_max = 2.5
            rKron_bins = np.arange(rKron_min,rKron_max,(rKron_max-rKron_min)/n_rKron_bins)
            mu,sigma = norm.fit(rKron)
            mu_pix = mu / wcs_pixscale_x
            x1,bins1,p1 = axPlot1.hist(rKron,bins=rKron_bins,histtype='bar',density=1,edgecolor='black',color='#2077b4',zorder=1)
            #y = mlab.normpdf(rKron_bins,mu,sigma)
            y = scipy.stats.norm.pdf(rKron_bins,mu,sigma)
            axPlot1.plot(rKron_bins,y,ls='--',lw=2,color='#880000',zorder=2)
            plt.savefig(filename + '.pdf',format='pdf',transparent=True)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'rKron histogram generation and Gaussian fitting (mu={:.2f} arcsec; sigma={:.2f}) completed successfully.'.format(mu,sigma))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: generate_rKron_histogram()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - FWHM histogram generation using generate_rKron_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - FWHM histogram generation using generate_rKron_histogram() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: generate_rKron_histogram()')
    return mu,mu_pix,sigma,keep_going


##### FUNCTION DEFINITIONS -- DATA EXPORTING #####

def write_photometry_data_toingest(base_path,data_path,photometry_data_toingest_filename,fits_filestem,mosaic_elem_table,headerinfo,sky_mean,sky_stddev,psf_width_mean_arcsec,source_density,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,limit_mag_ps,limit_mag_sb,thread_idx,keep_going,path_logfile,path_errorfile):
    # writes table of exposure data for later ingestion into database
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Writing exposure data to {:s} for ingestion...'.format(photometry_data_toingest_filename))
            with open(photometry_data_toingest_filename,'a') as f:
                filename_base = fits_filestem[:-4]
                ext_number    = int(fits_filestem[-3:])
                mosaic_element_id = -1
                for idx in range(0,len(mosaic_elem_table)):
                    if mosaic_elem_table[idx]['elem_idx'] == ext_number:
                        mosaic_element_id = mosaic_elem_table[idx]['element_id']
                filename_parts      = filename_base.split('-')
                base_filename       = filename_base+'.fits'
                if os.path.isfile(fits_filestem+'.photsolved'):
                    source_table_path   = data_path+'source_tables/'
                    source_table_file   = fits_filestem + '_calibrated_sources.txt'
                else:
                    source_table_path   = 'none'
                    source_table_file   = 'none'
                date_tai            = headerinfo[0]
                exposure_start_tai  = headerinfo[2]
                exposure_time       = float(headerinfo[3])
                filter_name         = headerinfo[4]
                pointing_center_ra  = float(headerinfo[6])
                pointing_center_dec = float(headerinfo[7])
                f.write('{:>17d}   {:<27s}   {:<10s}         {:s}  {:<11s}   {:15.10f}   {:16.10f}   {:16.10f}   {:15.10f}   {:21.3f}   {:14.3f}   {:14.3f}   {:14.3f}   {:>15d}   {:12.3f}   {:12.3f}   {:s}   {:s}\n'.format(mosaic_element_id,base_filename,date_tai,exposure_start_tai,filter_name,pointing_center_ra,pointing_center_dec,sky_mean,sky_stddev,psf_width_mean_arcsec,source_density,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,limit_mag_ps,limit_mag_sb,source_table_path,source_table_file))
            mcd.output_log_entry(path_logfile,'Writing exposure data done.')
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: write_photometry_data_toingest()'.format(fits_filestem))
            print(e)
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: write_photometry_data_toingest()')
    return keep_going


def write_calibrated_source_file(field_source_table,fits_filestem,avg_zeropoint_mag,avg_zeropoint_err,exptime,filter_name,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            num_sources = len(field_source_table)
            calibrated_source_filename = fits_filestem + '_calibrated_sources.txt'
            mcd.output_log_entry(path_logfile,'Writing calibrated sources data file...'.format(calibrated_source_filename))
            mcd.output_log_entry(path_logfile,'Number of sources: {:d}'.format(num_sources))
            with open(calibrated_source_filename,'w') as calibrated_source_file:
                mcd.output_log_entry(path_logfile,'Writing calibrated source file: {:s}'.format(calibrated_source_filename))
                calibrated_source_file.write('     x_s       y_s    flux_s   dflux_s  major_s  minor_s   phi_s    chiN_s      SNR_s       x_t       y_t    flux_t   dflux_t  trail_t   fwhm_t   phi_t   chiN_t      SNR_t          RA         Dec   filter       mag_s  magerr_s   mag_t  magerr_t\n')
                for idx in range(0,num_sources):
                    x_s     = field_source_table[idx]['x_s']
                    y_s     = field_source_table[idx]['y_s']
                    flux_s  = field_source_table[idx]['flux_s']
                    dflux_s = field_source_table[idx]['dflux_s']
                    major_s = field_source_table[idx]['major_s']
                    minor_s = field_source_table[idx]['minor_s']
                    phi_s   = field_source_table[idx]['phi_s']
                    chiN_s  = field_source_table[idx]['chiN_s']
                    SNR_s   = field_source_table[idx]['SNR_s']
                    x_t     = field_source_table[idx]['x_t']
                    y_t     = field_source_table[idx]['y_t']
                    flux_t  = field_source_table[idx]['flux_t']
                    dflux_t = field_source_table[idx]['dflux_t']
                    trail_t = field_source_table[idx]['trail_t']
                    fwhm_t  = field_source_table[idx]['fwhm_t']
                    phi_t   = field_source_table[idx]['phi_t']
                    chiN_t  = field_source_table[idx]['chiN_t']
                    SNR_t   = field_source_table[idx]['SNR_t']                    
                    RA      = field_source_table[idx]['RA']
                    Dec     = field_source_table[idx]['Dec']
                    mag_s,magerr_s,proc_status,keep_going = compute_mags_from_zpoint(flux_s,dflux_s,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                    mag_t,magerr_t,proc_status,keep_going = compute_mags_from_zpoint(flux_t,dflux_t,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                    calibrated_source_file.write('{:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}   {:6.2f}   {:6.2f}  {:6.2f}   {:7.2f}  {:9.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}   {:6.2f}   {:6.2f}  {:6.2f}  {:7.2f}  {:9.2f}  {:10.6f}  {:10.6f}   {:<10s}  {:6.3f}    {:6.3f}  {:6.3f}    {:6.3f}\n'.format(x_s,y_s,flux_s,dflux_s,major_s,minor_s,phi_s,chiN_s,SNR_s,x_t,y_t,flux_t,dflux_t,trail_t,fwhm_t,phi_t,chiN_t,SNR_t,RA,Dec,filter_name,mag_s,magerr_s,mag_t,magerr_t))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: write_calibrated_source_file()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Creation of calibrated source file using create_calibrated_source_file() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Creation of calibrated source file using create_calibrated_source_file() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Creation of calibrated source file using create_calibrated_source_file() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: write_calibrated_source_file()')
    return proc_status,keep_going


def write_multiap_photom_file(multiap_photom_table,fits_filestem,avg_zeropoint_mag,avg_zeropoint_err,exptime,filter_name,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            num_sources = len(multiap_photom_table)
            multiap_photom_filename = fits_filestem + '_multiap_photometry.txt'
            mcd.output_log_entry(path_logfile,'Writing multi-aperture photometry data file...'.format(multiap_photom_filename))
            mcd.output_log_entry(path_logfile,'Number of sources: {:d}'.format(num_sources))
            with open(multiap_photom_filename,'w') as multiap_photom_file:
                mcd.output_log_entry(path_logfile,'Writing multi-aperture photometry file: {:s}'.format(multiap_photom_filename))
                multiap_photom_file.write('     x_t       y_t     x_05rK    y_05rK  flux_t_05rK  dflux_t_05rK  mag_t_05rK  magerr_t_05rK     x_10rK    y_10rK  flux_t_10rK  dflux_t_10rK  mag_t_10rK  magerr_t_10rK     x_20rK    y_20rK  flux_t_20rK  dflux_t_20rK  mag_t_20rK  magerr_t_20rK     x_30rK    y_30rK  flux_t_30rK  dflux_t_30rK  mag_t_30rK  magerr_t_30rK     x_40rK    y_40rK  flux_t_40rK  dflux_t_40rK  mag_t_40rK  magerr_t_40rK     x_50rK    y_50rK  flux_t_50rK  dflux_t_50rK  mag_t_50rK  magerr_t_50rK\n')
                for idx in range(0,num_sources):
                    mag_t_05rK,magerr_t_05rK,mag_t_10rK,magerr_t_10rK,mag_t_20rK,magerr_t_20rK = 0,0,0,0,0,0
                    mag_t_30rK,magerr_t_30rK,mag_t_40rK,magerr_t_40rK,mag_t_50rK,magerr_t_50rK = 0,0,0,0,0,0
                    x_t          = multiap_photom_table[idx]['x_t']
                    y_t          = multiap_photom_table[idx]['y_t']
                    x_t_05rK     = multiap_photom_table[idx]['x_05rK']
                    y_t_05rK     = multiap_photom_table[idx]['y_05rK']
                    flux_t_05rK  = multiap_photom_table[idx]['flux_t_05rK']
                    dflux_t_05rK = multiap_photom_table[idx]['dflux_t_05rK']
                    x_t_10rK     = multiap_photom_table[idx]['x_10rK']
                    y_t_10rK     = multiap_photom_table[idx]['y_10rK']
                    flux_t_10rK  = multiap_photom_table[idx]['flux_t_10rK']
                    dflux_t_10rK = multiap_photom_table[idx]['dflux_t_10rK']
                    x_t_20rK     = multiap_photom_table[idx]['x_20rK']
                    y_t_20rK     = multiap_photom_table[idx]['y_20rK']
                    flux_t_20rK  = multiap_photom_table[idx]['flux_t_20rK']
                    dflux_t_20rK = multiap_photom_table[idx]['dflux_t_20rK']
                    x_t_30rK     = multiap_photom_table[idx]['x_30rK']
                    y_t_30rK     = multiap_photom_table[idx]['y_30rK']
                    flux_t_30rK  = multiap_photom_table[idx]['flux_t_30rK']
                    dflux_t_30rK = multiap_photom_table[idx]['dflux_t_30rK']
                    x_t_40rK     = multiap_photom_table[idx]['x_40rK']
                    y_t_40rK     = multiap_photom_table[idx]['y_40rK']
                    flux_t_40rK  = multiap_photom_table[idx]['flux_t_40rK']
                    dflux_t_40rK = multiap_photom_table[idx]['dflux_t_40rK']
                    x_t_50rK     = multiap_photom_table[idx]['x_50rK']
                    y_t_50rK     = multiap_photom_table[idx]['y_50rK']
                    flux_t_50rK  = multiap_photom_table[idx]['flux_t_50rK']
                    dflux_t_50rK = multiap_photom_table[idx]['dflux_t_50rK']
                    if flux_t_05rK != 0 and dflux_t_05rK != 0:
                        mag_t_05rK,magerr_t_05rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_05rK,dflux_t_05rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                    if flux_t_10rK != 0 and dflux_t_10rK != 0:
                        mag_t_10rK,magerr_t_10rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_10rK,dflux_t_10rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                    if flux_t_20rK != 0 and dflux_t_20rK != 0:
                        mag_t_20rK,magerr_t_20rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_20rK,dflux_t_20rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                    if flux_t_30rK != 0 and dflux_t_30rK != 0:
                        mag_t_30rK,magerr_t_30rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_30rK,dflux_t_30rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                    if flux_t_40rK != 0 and dflux_t_40rK != 0:
                        mag_t_40rK,magerr_t_40rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_40rK,dflux_t_40rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                    if flux_t_50rK != 0 and dflux_t_50rK != 0:
                        mag_t_50rK,magerr_t_50rK,proc_status,keep_going = compute_mags_from_zpoint(flux_t_50rK,dflux_t_50rK,exptime,avg_zeropoint_mag,avg_zeropoint_err,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                        keep_going = True
                        
                    multiap_photom_file.write('{:8.2f}  {:8.2f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}   {:8.2f}  {:8.2f}     {:8.2f}      {:8.2f}      {:6.3f}         {:6.3f}\n'.format(x_t,y_t, \
                        x_t_05rK,y_t_05rK,flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK, \
                        x_t_10rK,y_t_10rK,flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK, \
                        x_t_20rK,y_t_20rK,flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK, \
                        x_t_30rK,y_t_30rK,flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK, \
                        x_t_40rK,y_t_40rK,flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK, \
                        x_t_50rK,y_t_50rK,flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: write_multiap_photom_file()'.format(fits_filename_wcs))
            print(e)
            if os.path.exists(fits_filename_wcs[:-9]+'.photfailed'):
                with open(fits_filename_wcs[:-9]+'.photfailed','a') as f:
                    f.write('{:s} - Creation of multi-aperture photometry file using write_multiap_photom_file() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                with open(fits_filename_wcs[:-9]+'.photfailed','w') as f:
                    f.write('{:s} - Creation of multi-aperture photometry file using write_multiap_photom_file() failed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            proc_status = 'Creation of multi-aperture photometry file using write_multiap_photom_file() failed'
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: write_multiap_photom_file()')
    return proc_status,keep_going


##### FUNCTION DEFINITIONS -- CLEAN UP INTERMEDIATE FILES #####

def clean_phot_output_file(output_filename,sub_dir_path):
    if os.path.exists(output_filename):                # check if original file exists
        if os.path.exists(output_filename+'.gz'):      # remove previous compressed version of file in current directory if exists
            os.remove(output_filename+'.gz')
        if os.path.exists(sub_dir_path+output_filename):  # remove previous uncompressed version of file in sub-directory if exists
            os.remove(sub_dir_path+output_filename)
        if os.path.exists(sub_dir_path+output_filename+'.gz'):  # remove previous compressed version of file in sub-directory if exists
            os.remove(sub_dir_path+output_filename+'.gz')
        mcd.compress_file_gzip(output_filename)                     # compress original file
        os.rename(output_filename+'.gz',sub_dir_path+output_filename+'.gz')  # move compressed original file to sub-directory
    return None

def clean_phot_output(path_datadownload,fits_filestem,path_logfile,path_errorfile):
    try:
        mcd.output_log_entry(path_logfile,'Cleaning up photometric calibration output...')
        phot_output_path   = path_datadownload+'phot_output/'
        source_tables_path = path_datadownload+'source_tables/'
        mcd.create_directory(phot_output_path,path_logfile,path_errorfile)
        mcd.create_directory(source_tables_path,path_logfile,path_errorfile)
        
        clean_phot_output_file(fits_filestem+'_fwhm_histogram1.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_fwhm_histogram2.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_rKron_histogram1.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_rKron_histogram2.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_refcat_sources.dat',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.stars',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.moments',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.srclist',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.05rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.10rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.20rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.30rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.40rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_wcs.50rKron.trails',phot_output_path)
        clean_phot_output_file(fits_filestem+'_zeropoint_histogram1.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_zeropoint_histogram2.pdf',phot_output_path)
        clean_phot_output_file(fits_filestem+'_calibrated_sources.txt',source_tables_path)
        clean_phot_output_file(fits_filestem+'_multiap_photometry.txt',source_tables_path)        
        mcd.output_log_entry(path_logfile,'Cleaning up photometric calibration output done.')
        keep_going = True
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: clean_phot_output()'.format(fits_filestem))
    return None


def main():

    # Define filenames and paths
    if len(sys.argv)!=5:
        print('Usage:\n python3 pyt_astrometric_photometric_calibration_sdss.py [base_path] [sqlite_file] [instrument] [thread_idx]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    instrument  = sys.argv[3]
    thread_idx  = int(sys.argv[4])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('IM04b_photometric_calibration_file_list_{:s}_{:02d} execution started.'.format(instrument,thread_idx),'{:s} - IM04b_photometric_calibration_file_list_{:s}_{:02d} execution started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),instrument,thread_idx))

        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'IM04b_photometric_calibration_file_list_{:02d}'.format(thread_idx))
        
        os.chdir(base_path)
        dir_exposure_data = base_path+'exposure_data/'
        dir_files_toproc  = base_path+'ssois_queries/files_toproc/'
        mcd.create_directory(dir_exposure_data,path_logfile,path_errorfile)
        if instrument == 'SDSS':
            dir_processed_data = base_path+'data_sdss_processed/'
        #elif instrument == 'MegaCam':
        #    dir_processed_data = base_path+'data_megacam_processed/'
        else:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Instrument {:s} not recognized or not yet implemented.'.format(instrument))
            keep_going = False
            
        if not os.path.isdir(dir_processed_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Processed data directory {:s} not found.'.format(dir_processed_data))
            mcd.send_status_email('IM04b_photometric_calibration_file_list_{:02d} execution failed'.format(thread_idx),'IM04b_photometric_calibration_sdss_file_list_{:02d} execution for {:s} data failed - Processed data directory {:s} not found.'.format(thread_idx,instrument,dir_processed_data))
            keep_going = False
        if not os.path.isdir(dir_files_toproc):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'File lists directory {:s} not found.'.format(dir_processed_data))
            mcd.send_status_email('IM04b_photometric_calibration_file_list_{:02d} execution failed'.format(thread_idx),'IM04b_photometric_calibration_sdss_file_list_{:02d} execution for {:s} data failed - File lists directory {:s} not found.'.format(thread_idx,instrument,dir_processed_data))
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
            #instrument_id,keep_going = retrieve_instrument_id(sqlite_file,instrument_name,thread_idx,keep_going,path_logfile,path_errorfile)
            instrument_id = 1  # SDSS Imaging Camera
            mosaic_elem_table    = Table(rows=mosaic_elem_datarows,names=('elem_idx','element_id'))
            for idx in range(0,num_elements):
                mosaic_element_num = mosaic_elem_table[idx]['elem_idx']
                #mosaic_element_id,keep_going = retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile)
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
            # initialize exposure data output file
            photometry_data_toingest_filename = dir_exposure_data+'00_photometry_data_{:s}_{:s}_{:02d}_toingest.txt'.format(instrument,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),thread_idx)
            processed_files_toingest_filename = dir_exposure_data+'00_processed_files_photometry_{:s}_{:s}_{:02d}_toingest.txt'.format(instrument,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),thread_idx)
            mcd.output_log_entry(path_logfile,'Writing photometric output to {:s}...'.format(photometry_data_toingest_filename))
            mcd.output_log_entry(path_logfile,'Writing processing statuses to {:s}...'.format(processed_files_toingest_filename))
            if not os.path.exists(photometry_data_toingest_filename):
                with open(photometry_data_toingest_filename,'w') as f:
                    f.write('mosaic_element_id   base_filename                   date_tai   exposure_start_tai  filter_name   pointing_ctr_ra   pointing_ctr_dec           sky_mean        sky_stddev   psf_width_mean_arcsec   source_density   avg_zpoint_mag   avg_zpoint_err   num_calib_stars   limit_mag_ps   limit_mag_sb   source_table_path   source_table_file\n')
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
                        for fits_filename_fz in sorted(glob.glob(filename_toprocess[:-9]+'_*.fits.fz')):
                            if os.path.isfile(fits_filename_fz[:-12]+'.photsolved') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.photfailed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.sdss_header_extraction_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.radec_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.lowsnr_source_removal_stopped') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.lowsnr_source_removal_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.skylevel_measurement_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.source_extraction_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.field_source_table_creation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.bad_source_removal_stopped') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.bad_source_removal_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.lowsnr_source_removal_stopped') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.lowsnr_source_removal_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.source_density_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.detection_limits_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.refcat_source_retrieval_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.refcat_source_matching_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.detection_limits_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.multiap_photometry_tphot_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.kron_radius_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.psf_width_computation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.calibrated_source_file_creation_failed') \
                                or os.path.isfile(fits_filename_fz[:-12]+'.multiap_photometry_file_creation_failed'):
                                mcd.output_log_entry(path_logfile,'{:s} already processed.'.format(fits_filename_fz))
                                phot_done = True
                            else:
                                phot_done = False
                                if fits_filename_fz[-12:-8] == '_wcs' and fits_filename_fz[6:7] != 'u':
                                    keep_going = True
                                    mcd.output_log_entry(path_logfile,'Decompressing {:s} for processing...'.format(fits_filename_fz))
                                    proc_status = ''
                                    mcd.decompress_file_funpack(fits_filename_fz)
                                    fits_filename = fits_filename_fz[:-3]  # e.g., frame-z-004828-1-0430_000_wcs.fits.fz --> frame-z-004828-1-0430_000_wcs.fits
                                    fits_filestem = fits_filename[:-9]  # e.g., frame-z-004828-1-0430_000_wcs.fits --> frame-z-004828-1-0430_000
                                    fits_filename_wcs = fits_filestem + '_wcs.fits'  # e.g., frame-z-004828-1-0430_000 --> frame-z-004828-1-0430_000_wcs.fits
                                        
                                    ### Extract header info from FITS file
                                    if instrument == 'SDSS':
                                        headerinfo,proc_status,keep_going = extract_header_info_sdss(fits_filename,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    exptime     = float(headerinfo[3])
                                    filter_name = headerinfo[4]
                                    exp_npix_x  = headerinfo[12]
                                    exp_npix_y  = headerinfo[13]
                                    ra_deg      = headerinfo[14]
                                    dec_deg     = headerinfo[15]
                                    wcs_pixscale_x,wcs_pixscale_y = compute_wcs_pixel_scale(fits_filename)
                                    radius_arcmin  = exp_npix_x*wcs_pixscale_x/60

                                    ### Compute photometric calibration solution for FITS file
                                    sky_mean,sky_stddev,proc_status,keep_going = measure_skylevel_stddev(fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)

                                    # Source extraction
                                    # 1. extract_sources_sdss():              Source extraction for target identification and photometry (using lower S/N threshold)
                                    # 2. create_field_source_table():         Create and populate field source table
                                    # 3. remove_bad_sources():                Remove bad sources (where flux < 0 or SNR < 3)
                                    # 4. compute_exp_source_density():        Compute overall source density of exposure
                                    # 5. compute_psf_width_mean():            Compute mean PSF width of field sources with outlier rejection
                                    # 6. get_radec_from_pixel_coords_array(): Convert pixel coordinates of extracted catalog sources to RA/Dec

                                    if instrument == 'SDSS':
                                        output_path_stars,output_path_trails,output_path_moments,output_path_source_list,proc_status,keep_going = extract_sources_sdss(sky_stddev,fits_filename_wcs,fits_filestem,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    field_source_table,calib_source_table,proc_status,keep_going = create_field_source_table(output_path_stars,output_path_trails,output_path_moments,fits_filename_wcs,match_threshold_pix,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    field_source_table,proc_status,keep_going = remove_bad_sources(field_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    source_density,proc_status,keep_going     = compute_exp_source_density(exp_npix_x,exp_npix_y,wcs_pixscale_x,wcs_pixscale_y,field_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    psf_width_mean_arcsec,psf_width_mean_pix,calib_source_table,proc_status,keep_going = compute_psf_width_mean(wcs_pixscale_x,calib_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    rKron_mean_arcsec,rKron_mean_pix,calib_source_table,proc_status,keep_going = compute_rKron_mean(wcs_pixscale_x,calib_source_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    field_source_table,proc_status,keep_going = get_radec_from_pixel_coords_array(fits_filename_wcs,field_source_table,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    calib_source_table,proc_status,keep_going = get_radec_from_pixel_coords_array(fits_filename_wcs,calib_source_table,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)

                                    # Perform photometric calibration
                                    # 1. retrieve_refcat_sources():  Retrieve refcat sources for photometric calibration
                                    # 2. match_refcat_sources():     Match extracted field sources with refcat sources
                                    # 3. remove_lowsnr_sources():    Remove sources with SNR < 3
                                    # 4. compute_zero_point():       Compute average zero-point of image
                                    # 5. compute_detection_limits(): Compute point-source and surface brightness detection limits of image
                                    # 6. clean_phot_output():        Clean up photometric calibration files
        
                                    refcat_table,proc_status,keep_going = retrieve_refcat_sources(ra_deg,dec_deg,radius_arcmin,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    calib_star_table,proc_status,keep_going = match_refcat_sources(calib_source_table,refcat_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    calib_star_table,proc_status,keep_going = remove_lowsnr_sources(calib_star_table,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    calib_star_table,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,proc_status,keep_going = compute_zero_point(fits_filename_wcs,filter_name,calib_star_table,exptime,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    limit_mag_ps,limit_mag_sb,proc_status,keep_going = compute_detection_limits(sky_stddev,psf_width_mean_pix,wcs_pixscale_x,exptime,avg_zeropoint_mag,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    proc_status,keep_going = write_calibrated_source_file(field_source_table,fits_filestem,avg_zeropoint_mag,avg_zeropoint_err,exptime,filter_name,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    multiap_photom_table,proc_status,keep_going = perform_multiap_photometry_tphot(instrument,field_source_table,rKron_mean_pix,sky_stddev,output_path_source_list,fits_filename_wcs,fits_filestem,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    proc_status,keep_going = write_multiap_photom_file(multiap_photom_table,fits_filestem,avg_zeropoint_mag,avg_zeropoint_err,exptime,filter_name,fits_filename_wcs,proc_status,thread_idx,keep_going,path_logfile,path_errorfile)
                                    keep_going = write_photometry_data_toingest(base_path,data_path,photometry_data_toingest_filename,fits_filestem,mosaic_elem_table,headerinfo,sky_mean,sky_stddev,psf_width_mean_arcsec,source_density,avg_zeropoint_mag,avg_zeropoint_err,num_calib_stars,limit_mag_ps,limit_mag_sb,thread_idx,keep_going,path_logfile,path_errorfile)
                                    clean_phot_output(data_path,fits_filestem,path_logfile,path_errorfile)

                                    with open(processed_files_toingest_filename,'a') as of:
                                        of.write('{:<50s}   {:s}\n'.format(fits_filename_fz,proc_status))
                                        
                                    # Compress fits file for archiving
                                    mcd.output_log_entry(path_logfile,'Re-compressing {:s}...'.format(fits_filename_wcs))
                                    mcd.compress_file_fpack(fits_filename_wcs)

                                else:
                                    mcd.output_log_entry(path_logfile,'Astrometric calibration for {:s} not yet complete.'.format(fits_filename_fz))
                                    keep_going = False
                    else:
                        mcd.output_error_log_entry(path_logfile,path_errorfile,'Data path {:s} not found.'.format(data_path))
                        
            mcd.compress_file_gzip(photometry_data_toingest_filename)
            mcd.compress_file_gzip(processed_files_toingest_filename)
            mcd.compress_file_gzip(ssois_query_results_filepath)

        mcd.output_log_entry(path_logfile,'IM04b_photometric_calibration_file_list_{:s}_{:02d} execution complete.'.format(instrument,thread_idx))
        mcd.send_status_email('IM04b_photometric_calibration_file_list_{:s}_{:02d} execution complete.'.format(instrument,thread_idx),'IM04b_photometric_calibration_file_list_{:s}_{:02d} execution complete.'.format(instrument,thread_idx))

    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None


if __name__ == '__main__':
    main()

