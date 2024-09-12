import os
import os.path
import sys
import re
import glob
import subprocess
import datetime
import math
import ssl
import bz2

import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
import numpy.ma as ma
import scipy
import statistics
import uncertainties
import urllib

from astropy import units as u
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord,Angle
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.stats import SigmaClip,sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import simple_norm
from astroquery.jplhorizons import Horizons
from bz2 import decompress
from matplotlib.ticker import NullFormatter, MaxNLocator, ScalarFormatter, MultipleLocator, FormatStrFormatter
from photutils.background import Background2D, MedianBackground
from photutils.aperture import CircularAperture,CircularAnnulus,aperture_photometry
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry
from scipy.stats import norm
from uncertainties import unumpy,ufloat
from uncertainties.umath import *
from uncertainties.unumpy import *


########## TOP LEVEL FUNCTION ##########

def search_extract_archival_photometry(param_dict):
    os.chdir(param_dict['base_path'])

    # Perform setup for processing (making data directories, loading config parameters, etc.)
    param_dict = initialize_processing(param_dict)

    # Run search for potential detections in archival data using SSOIS
    param_dict = run_ssois_search(param_dict)
    
    # Download and perform initial processing of image files identified by SSOIS search
    param_dict = download_ssois_query_result_files(param_dict)
    
    # Locate objects in images and perform multi-aperture photometry
    param_dict = perform_photometry(param_dict)

    # Move final results files into target photometry folder
    param_dict = clean_final_photometry_files(param_dict)
    
    return None


########## INITIAL PROCESSING FUNCTIONS ##########

def initialize_processing(param_dict):
    param_dict = set_default_photometry_params(param_dict)
    param_dict = set_epochs(param_dict)
    param_dict = create_ssois_search_directory(param_dict)
    param_dict['log_filepath'] = initialize_log_file(param_dict['query_results_dirpath'])
    param_dict = write_phot_params_to_log(param_dict)
    param_dict = read_config_file(param_dict)
    param_dict = get_instrument_params(param_dict)
    return param_dict

def set_default_photometry_params(param_dict):
    if 'field_source_phot_radius_arcsec' not in param_dict:
        param_dict['field_source_phot_radius_arcsec'] = 3.5
    if 'sky_inner_r_arcsec' not in param_dict:
        param_dict['sky_inner_r_arcsec'] = 15
    if 'sky_outer_r_arcsec' not in param_dict:
        param_dict['sky_outer_r_arcsec'] = 25
    if 'daofind_fwhm' not in param_dict:
        param_dict['daofind_fwhm'] = 1.5
    if 'daofind_sigma_threshold' not in param_dict:
        param_dict['daofind_sigma_threshold'] = 5
    if 'daofind_fwhm_target' not in param_dict:
        param_dict['daofind_fwhm_target'] = 1.5
    if 'daofind_sigma_threshold_target' not in param_dict:
        param_dict['daofind_sigma_threshold_target'] = 3
    if 'target_sky_inner_r_arcsec' not in param_dict:
        param_dict['target_sky_inner_r_arcsec'] = 15
    if 'target_sky_outer_r_arcsec' not in param_dict:
        param_dict['target_sky_outer_r_arcsec'] = 25
    if 'target_phot_radii_arcsec' not in param_dict:
        param_dict['target_phot_radii_arcsec'] = np.arange(1,10.5,0.5)
    return param_dict

def set_epochs(param_dict):
    param_dict['epoch1'] = '1990-01-01'
    param_dict['epoch2'] = datetime.datetime.today().strftime('%Y-%m-%d')
    if 'start_date' in param_dict:
        if param_dict['start_date'] != '':
            param_dict['epoch1'] = param_dict['start_date']
    if 'end_date' in param_dict:
        if param_dict['end_date'] != '':
            param_dict['epoch2'] = param_dict['end_date']
    return param_dict

def create_ssois_search_directory(param_dict):
    target_name = remove_extra_characters(param_dict['target_name'])
    telinst     = remove_extra_characters(param_dict['telinst'])
    epoch1      = remove_extra_characters(param_dict['epoch1'])
    epoch2      = remove_extra_characters(param_dict['epoch2'])
    obj_results_dirpath   = param_dict['base_path'] + 'archival_search_{:s}/'.format(target_name)
    query_results_dirpath = obj_results_dirpath + '{:s}_{:s}_{:s}/'.format(telinst,epoch1,epoch2)
    if not os.path.exists(obj_results_dirpath):
        os.mkdir(obj_results_dirpath)
    if not os.path.exists(query_results_dirpath):
        os.mkdir(query_results_dirpath)
    param_dict['query_results_dirpath']  = query_results_dirpath
    return param_dict

def initialize_log_file(base_path):
    # Create and open log file
    log_filepath = base_path + 'log_archival_photometry_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    with open(log_filepath,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return log_filepath

def write_phot_params_to_log(param_dict):
    output_log_entry(param_dict['log_filepath'],'target:          {:s}'.format(param_dict['target_name']))
    output_log_entry(param_dict['log_filepath'],'telescope:       {:s}'.format(param_dict['telinst']))
    output_log_entry(param_dict['log_filepath'],'field_source_phot_radius_arcsec: {:.2f}'.format(param_dict['field_source_phot_radius_arcsec']))
    output_log_entry(param_dict['log_filepath'],'sky_inner_r_arcsec:              {:.2f}'.format(param_dict['sky_inner_r_arcsec']))
    output_log_entry(param_dict['log_filepath'],'sky_outer_r_arcsec:              {:.2f}'.format(param_dict['sky_outer_r_arcsec']))
    output_log_entry(param_dict['log_filepath'],'daofind_fwhm:                    {:.2f}'.format(param_dict['daofind_fwhm']))
    output_log_entry(param_dict['log_filepath'],'daofind_fwhm_target:             {:.2f}'.format(param_dict['daofind_fwhm_target']))
    output_log_entry(param_dict['log_filepath'],'daofind_sigma_threshold:         {:.2f}'.format(param_dict['daofind_sigma_threshold']))
    for idx in range(len(param_dict['target_phot_radii_arcsec'])):
        output_log_entry(param_dict['log_filepath'],'target_phot_radius_arcsec:       {:.2f}'.format(param_dict['target_phot_radii_arcsec'][idx]))
    output_log_entry(param_dict['log_filepath'],'target_sky_inner_r_arcsec:       {:.2f}'.format(param_dict['target_sky_inner_r_arcsec']))
    output_log_entry(param_dict['log_filepath'],'target_sky_outer_r_arcsec:       {:.2f}'.format(param_dict['target_sky_outer_r_arcsec']))
    return param_dict

def read_config_file(param_dict):
    config = {}
    config_filepath = param_dict['base_path'] + 'macadamia.cfg'
    with open(config_filepath,'r') as config_file:
        for line in config_file:
            line_data = line.split(',')
            config[line_data[0]] = line_data[1].strip()
    param_dict['config'] = config
    return param_dict

def get_instrument_params(param_dict):
    instrument_found = False
    instrument_params = {}
    instrument = param_dict['telinst']
    if instrument == 'SDSS':
        instrument_params['pixscale']      = 0.396
        instrument_params['first_element'] = 0
        instrument_params['num_elements']  = 1
        instrument_params['obs_code']      = '645'
        instrument_found = True
    if instrument == 'CTIO-4m/DECam':
        instrument_params['pixscale']      = 0.2637
        instrument_params['first_element'] = 1
        instrument_params['num_elements']  = 60
        instrument_params['obs_code']      = '807'
        instrument_found = True
    if instrument == 'CFHT/MegaCam':
        instrument_params['pixscale']      = 0.187
        instrument_params['first_element'] = 1
        instrument_params['num_elements']  = 40
        instrument_params['obs_code']      = '568'
        instrument_found = True
    if instrument_found:
        output_log_entry(param_dict['log_filepath'],'{:s} found:'.format(instrument))
        output_log_entry(param_dict['log_filepath'],' - Pixel scale: {:f}'.format(instrument_params['pixscale']))
        output_log_entry(param_dict['log_filepath'],' - Index of first element: {:d}'.format(instrument_params['first_element']))
        output_log_entry(param_dict['log_filepath'],' - Number of elements: {:d}'.format(instrument_params['num_elements']))
    else:
        output_log_entry(param_dict['log_filepath'],'{:s} not found:'.format(instrument))
    param_dict['instrument_params'] = instrument_params
    return param_dict


########## SSOIS SEARCH FUNCTIONS ##########

def run_ssois_search(param_dict):
    # Run search for potential detections in archival data using SSOIS
    param_dict = generate_ssois_search_url_results_file(param_dict)
    param_dict = execute_ssois_search(param_dict)
    if param_dict['telinst'] == 'CTIO-4m/DECam':
        param_dict = parse_ssois_query_results_decam(param_dict)
    else:
        param_dict = parse_ssois_query_results(param_dict)
    param_dict = run_horizons_queries_prelim(param_dict)
    return param_dict

def generate_ssois_search_url_results_file(param_dict):
    target_name = param_dict['target_name'].replace(' ','+')
    target_name = target_name.replace('/','%2F')
    telinst     = param_dict['telinst']
    epoch1      = param_dict['epoch1']
    epoch2      = param_dict['epoch2']
    param_dict['ssois_search_url'] = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl?lang=en;object={:s};search=bynameall;epoch1={:s};epoch2={:s};eellipse=;eunits=arcseconds;extres=no;xyres=no;telinst={:s};format=tsv'.format(target_name,epoch1,epoch2,telinst)
    output_log_entry(param_dict['log_filepath'],'SSOIS search URL: {:s}'.format(param_dict['ssois_search_url']))

    query_results_filename = 'ssois_results_{:s}_{:s}_{:s}_{:s}.txt'.format(remove_extra_characters(param_dict['target_name']),
        remove_extra_characters(param_dict['telinst']),remove_extra_characters(param_dict['epoch1']),remove_extra_characters(param_dict['epoch2']))
    param_dict['query_results_filepath'] = param_dict['query_results_dirpath'] + query_results_filename
    output_log_entry(param_dict['log_filepath'],'SSOIS query results filepath: {:s}'.format(param_dict['query_results_filepath']))

    return param_dict

def execute_ssois_search(param_dict):
    url = param_dict['ssois_search_url']
    result_filepath = param_dict['query_results_filepath']
    if not os.path.exists(result_filepath):
        output_log_entry(param_dict['log_filepath'],'Starting SSOIS query...')
        query_successful = True
        max_retries = 50
        for _ in range(max_retries):
            try:
                urllib.request.urlretrieve(url.strip(),result_filepath)
            except Exception:
                output_log_entry(param_dict['log_filepath'],'SSOIS query failed. Retrying...'.format(url))
            else:
                break
        else:
            output_log_entry(param_dict['log_filepath'],'SSOIS query failed ({:s}). Maximum retries reached'.format(url))
            query_successful = False
        if query_successful:
            output_log_entry(param_dict['log_filepath'],'Saving SSOIS query results to {:s}'.format(result_filepath))
    else:
        output_log_entry(param_dict['log_filepath'],'SSOIS query results file already found.')
    return param_dict

def parse_ssois_query_results(param_dict):
    query_results_filepath = param_dict['query_results_filepath']
    query_results = {}
    with open(query_results_filepath,'r') as query_results_file:
        for _ in range(1): #skip first header line
            next(query_results_file)
        for line in query_results_file:
            line_data  = re.split(r'\t+',line.rstrip('\t'))
            image_name = line_data[0]
            query_results[image_name] = {}
            query_results[image_name]['mjd']         = float(str(line_data[1]))
            query_results[image_name]['filter_name'] = line_data[2]
            query_results[image_name]['exptime']     = float(str(line_data[3]))
            query_results[image_name]['obj_ra']      = float(str(line_data[4]))
            query_results[image_name]['obj_dec']     = float(str(line_data[5]))
            query_results[image_name]['img_target']  = line_data[6]
            query_results[image_name]['datalink']    = line_data[8].strip()
    param_dict['query_results'] = query_results
    output_log_entry(param_dict['log_filepath'],'Query results read in from {:s}'.format(query_results_filepath))
    return param_dict

def parse_ssois_query_results_decam(param_dict):
    query_results_filepath = param_dict['query_results_filepath']
    query_results = {}
    prev_image_name,prev_image_prefix,prev_image_suffix = '','',''
    with open(query_results_filepath,'r') as query_results_file:
        for _ in range(1): #skip first header line
            next(query_results_file)
        for line in query_results_file:
            line_data  = re.split(r'\t+',line.rstrip('\t'))
            image_name = line_data[0].split()[0].strip()
            if '_ooi_' in image_name:
                image_prefix = image_name[:17]
                image_suffix = image_name[24:-8]
                if image_prefix == prev_image_prefix:
                    if prev_image_suffix != 'v1':
                        if image_suffix == 'v1' or ('v' in image_suffix and 'v' not in prev_image_suffix) or ('ls' in image_suffix and 'a' in prev_image_suffix):
                            del query_results[prev_image_name]
                            query_results[image_name] = {}
                            query_results[image_name]['mjd']         = line_data[1]
                            query_results[image_name]['filter_name'] = line_data[2]
                            query_results[image_name]['exptime']     = line_data[3]
                            query_results[image_name]['obj_ra']      = line_data[4]
                            query_results[image_name]['obj_dec']     = line_data[5]
                            query_results[image_name]['img_target']  = line_data[6]
                            query_results[image_name]['datalink']    = line_data[8].strip()
                            prev_image_name   = image_name
                            prev_image_prefix = image_prefix
                            prev_image_suffix = image_suffix
                else:
                    query_results[image_name] = {}
                    query_results[image_name]['mjd']         = line_data[1]
                    query_results[image_name]['filter_name'] = line_data[2]
                    query_results[image_name]['exptime']     = line_data[3]
                    query_results[image_name]['obj_ra']      = line_data[4]
                    query_results[image_name]['obj_dec']     = line_data[5]
                    query_results[image_name]['img_target']  = line_data[6]
                    query_results[image_name]['datalink']    = line_data[8].strip()
                    prev_image_name   = image_name
                    prev_image_prefix = image_prefix
                    prev_image_suffix = image_suffix
    param_dict['query_results'] = query_results
    output_log_entry(param_dict['log_filepath'],'Query results read in from {:s}'.format(query_results_filepath))
    return param_dict


########## DATA DOWNLOADING AND INITIAL PROCESSING FUNCTIONS ##########

def download_ssois_query_result_files(param_dict):
    output_log_entry(param_dict['log_filepath'],'Downloading {:s} data for {:s} from {:s} to {:s}...'.format(param_dict['telinst'],param_dict['target_name'],param_dict['epoch1'].replace('+','-'),param_dict['epoch2'].replace('+','-')))
    query_results = param_dict['query_results']
    images_to_delete = []
    for image_name in query_results:
        download_filename = ''
        if param_dict['telinst'] == 'CTIO-4m/DECam':
            if float(str(param_dict['query_results'][image_name]['exptime'])) >= 30 and param_dict['query_results'][image_name]['filter_name'] not in ['u','Y']:
                download_filename = param_dict['query_results_dirpath'] + image_name
            else:
                images_to_delete.append(image_name)
        elif param_dict['telinst'] == 'CFHT/MegaCam':
            if float(str(param_dict['query_results'][image_name]['exptime'])) >= 30 and param_dict['query_results'][image_name]['filter_name'] not in ['u','Y']:
                download_filename = param_dict['query_results_dirpath'] + os.path.basename(query_results[image_name]['datalink'])
            else:
                images_to_delete.append(image_name)
        else:
            download_filename = param_dict['query_results_dirpath'] + os.path.basename(query_results[image_name]['datalink'])

        if download_filename != '' and not os.path.exists(download_filename) and not os.path.exists(download_filename[:-3]) and not os.path.exists(download_filename[:-4]):
            try:
                output_log_entry(param_dict['log_filepath'],'Downloading {:s} to {:s}'.format(query_results[image_name]['datalink'],download_filename))
                urllib.request.urlretrieve(query_results[image_name]['datalink'],download_filename)
            except:
                output_log_entry(param_dict['log_filepath'],'{:s} could not be downloaded ({:s})'.format(download_filename,query_results[image_name]['datalink']))
                images_to_delete.append(image_name)
    for image_name in images_to_delete:
        del param_dict['query_results'][image_name]
    param_dict = decompress_data_get_metadata(param_dict)
    return param_dict

def decompress_data_get_metadata(param_dict):
    if param_dict['telinst'] == 'SDSS':
        param_dict = decompress_data_bz2(param_dict)
        param_dict = extract_metadata_sdss(param_dict)
    if param_dict['telinst'] == 'CTIO-4m/DECam':
        param_dict = decompress_data_fz(param_dict)
        param_dict = extract_metadata_decam(param_dict)
    if param_dict['telinst'] == 'CFHT/MegaCam':
        param_dict = decompress_data_fz(param_dict)
        param_dict = extract_metadata_megacam(param_dict)
    return param_dict

def extract_metadata_sdss(param_dict):
    for fits_filename in param_dict['fits_filenames']:
        with fits.open(fits_filename) as hdulist:
            hdr  = hdulist[0].header
            date_utc       = hdr['DATE-OBS']
            time_start_utc = hdr['TAIHMS']
            exposure_time  = float(hdr['EXPTIME'])
            time_start     = Time('{:s} {:s}'.format(date_utc,time_start_utc),scale='utc',format='iso')
            time_start_jd  = float(time_start.jd)
            mid_dt         = TimeDelta(exposure_time,format='sec') / 2
            time_mid       = time_start + mid_dt
            time_mid_jd    = float(time_mid.jd)
            filter_name    = hdr['FILTER'].strip()
            print('filter_name for {:s}: {:s} --> '.format(fits_filename,filter_name),end='')
            if filter_name == 'u': filter_name = 'u_sdss'
            if filter_name == 'g': filter_name = 'g_sdss'
            if filter_name == 'r': filter_name = 'r_sdss'
            if filter_name == 'i': filter_name = 'i_sdss'
            if filter_name == 'z': filter_name = 'z_sdss'
            print('{:s}'.format(filter_name))
            param_dict['fits_filenames'][fits_filename]['date_utc']       = date_utc
            param_dict['fits_filenames'][fits_filename]['time_start_utc'] = time_start_utc
            param_dict['fits_filenames'][fits_filename]['exposure_time']  = exposure_time
            param_dict['fits_filenames'][fits_filename]['time_start']     = time_start
            param_dict['fits_filenames'][fits_filename]['time_start_jd']  = time_start_jd
            param_dict['fits_filenames'][fits_filename]['time_mid']       = time_mid
            param_dict['fits_filenames'][fits_filename]['time_mid_jd']    = time_mid_jd
            param_dict['fits_filenames'][fits_filename]['filter']         = filter_name
    return param_dict

def extract_metadata_decam(param_dict):
    images_to_delete = []
    for fits_filename in param_dict['fits_filenames']:
        try:
            with fits.open(fits_filename) as hdulist:
                hdr  = hdulist[0].header
                date_utc       = hdr['DATE-OBS'][:10]
                time_start_utc = hdr['DATE-OBS'][11:]
                exposure_time  = float(hdr['EXPTIME'])
                time_start     = Time('{:s} {:s}'.format(date_utc,time_start_utc),scale='utc',format='iso')
                time_start_jd  = float(time_start.jd)
                mid_dt         = TimeDelta(exposure_time,format='sec') / 2
                time_mid       = time_start + mid_dt
                time_mid_jd    = float(time_mid.jd)
                filter_name    = hdr['FILTER'].strip()
                print('filter_name for {:s}: {:s} --> '.format(fits_filename,filter_name),end='')
                if filter_name == 'g DECam SDSS c0001 4720.0 1520.0': filter_name = 'g_sdss'
                if filter_name == 'r DECam SDSS c0002 6415.0 1480.0': filter_name = 'r_sdss'
                if filter_name == 'i DECam SDSS c0003 7835.0 1470.0': filter_name = 'i_sdss'
                if filter_name == 'z DECam SDSS c0004 9260.0 1520.0': filter_name = 'z_sdss'
                if filter_name == 'Y DECam c0005 10095.0 1130.0':     filter_name = 'Y'
                if filter_name == 'VR DECam c0007 6300.0 2600.0':     filter_name = 'VR'
                print('{:s}'.format(filter_name))
                param_dict['fits_filenames'][fits_filename]['date_utc']       = date_utc
                param_dict['fits_filenames'][fits_filename]['time_start_utc'] = time_start_utc
                param_dict['fits_filenames'][fits_filename]['exposure_time']  = exposure_time
                param_dict['fits_filenames'][fits_filename]['time_start']     = time_start
                param_dict['fits_filenames'][fits_filename]['time_start_jd']  = time_start_jd
                param_dict['fits_filenames'][fits_filename]['time_mid']       = time_mid
                param_dict['fits_filenames'][fits_filename]['time_mid_jd']    = time_mid_jd
                param_dict['fits_filenames'][fits_filename]['filter']         = filter_name
        except:
            output_log_entry(param_dict['log_filepath'],'Extraction of metadata from {:s} failed.'.format(fits_filename))
            images_to_delete.append(fits_filename)
    for fits_filename in images_to_delete:
        if fits_filename in param_dict['fits_filenames']:
            del param_dict['fits_filenames'][fits_filename]
    return param_dict

def extract_metadata_megacam(param_dict):
    images_to_delete = []
    for fits_filename in param_dict['fits_filenames']:
        try:
            with fits.open(fits_filename) as hdulist:
                hdr  = hdulist[0].header
                date_utc       = hdr['DATE-OBS']
                time_start_utc = hdr['UTC-OBS']
                exposure_time  = float(hdr['EXPTIME'])
                time_start     = Time('{:s} {:s}'.format(date_utc,time_start_utc),scale='utc',format='iso')
                time_start_jd  = float(time_start.jd)
                mid_dt         = TimeDelta(exposure_time,format='sec') / 2
                time_mid       = time_start + mid_dt
                time_mid_jd    = float(time_mid.jd)
                filter_name    = hdr['FILTER'].strip()
                print('filter_name for {:s}: {:s} --> '.format(fits_filename,filter_name),end='')
                if filter_name == 'u.MP9301':   filter_name = 'u_sdss'
                if filter_name == 'g.MP9401':   filter_name = 'g_sdss'
                if filter_name == 'g.MP9402':   filter_name = 'g_sdss'
                if filter_name == 'r.MP9601':   filter_name = 'r_sdss'
                if filter_name == 'r.MP9602':   filter_name = 'r_sdss'
                if filter_name == 'i.MP9701':   filter_name = 'i_sdss'
                if filter_name == 'i.MP9702':   filter_name = 'i_sdss'
                if filter_name == 'i.MP9703':   filter_name = 'i_sdss'
                if filter_name == 'z.MP9801':   filter_name = 'z_sdss'
                if filter_name == 'z.MP9901':   filter_name = 'z_sdss'
                if filter_name == 'gri.MP9605': filter_name = 'gri'
                print('{:s}'.format(filter_name))
                param_dict['fits_filenames'][fits_filename]['date_utc']       = date_utc
                param_dict['fits_filenames'][fits_filename]['time_start_utc'] = time_start_utc
                param_dict['fits_filenames'][fits_filename]['exposure_time']  = exposure_time
                param_dict['fits_filenames'][fits_filename]['time_start']     = time_start
                param_dict['fits_filenames'][fits_filename]['time_start_jd']  = time_start_jd
                param_dict['fits_filenames'][fits_filename]['time_mid']       = time_mid
                param_dict['fits_filenames'][fits_filename]['time_mid_jd']    = time_mid_jd
                param_dict['fits_filenames'][fits_filename]['filter']         = filter_name
        except:
            output_log_entry(param_dict['log_filepath'],'Extraction of metadata from {:s} failed.'.format(fits_filename))
            images_to_delete.append(fits_filename)
    for fits_filename in images_to_delete:
        if fits_filename in param_dict['fits_filenames']:
            del param_dict['fits_filenames'][fits_filename]
    return param_dict


########## COMPRESSION/DECOMPRESSION FUNCTIONS ##########

def decompress_data_bz2(param_dict):
    output_log_entry(param_dict['log_filepath'],'Decompressing .bz2 data...')
    os.chdir(param_dict['query_results_dirpath'])
    param_dict['fits_filenames'] = {}
    for image_name in param_dict['query_results']:
        filename_bz2 = os.path.basename(param_dict['query_results'][image_name]['datalink'])
        filename = filename_bz2[:-4]
        param_dict['fits_filenames'][filename] = {}
        if os.path.exists(filename_bz2) and not os.path.exists(filename):
            with bz2.open(filename_bz2,'rb') as f, open(filename,'wb') as of:
                data = f.read()
                of.write(data)
            os.remove(filename_bz2)
    output_log_entry(param_dict['log_filepath'],'Decompressing .bz2 data complete.')
    return param_dict

def decompress_data_fz(param_dict):
    output_log_entry(param_dict['log_filepath'],'Decompressing .fz data...')
    os.chdir(param_dict['query_results_dirpath'])
    param_dict['fits_filenames'] = {}
    images_to_delete = []
    for image_name in param_dict['query_results']:
        if param_dict['telinst'] == 'CTIO-4m/DECam':
            filename_fz = image_name
        else:
            filename_fz = os.path.basename(param_dict['query_results'][image_name]['datalink'])
        filename = filename_fz[:-3]
        param_dict['fits_filenames'][filename] = {}
        try:
            if os.path.exists(filename_fz) and not os.path.exists(filename):
                funpack(filename_fz,param_dict['config'])
        except:
            del param_dict['fits_filenames'][filename]
    output_log_entry(param_dict['log_filepath'],'Decompressing .fz data complete.')
    return param_dict


########## GENERAL PHOTOMETRY FUNCTIONS ##########

def perform_photometry(param_dict):
    filters = ['g_sdss','r_sdss','i_sdss','z_sdss','B','V','R','I']
    param_dict = convert_field_source_apert_params_to_pixels(param_dict)
    param_dict = convert_target_apert_params_to_pixels(param_dict)
    param_dict = run_horizons_queries(param_dict)
    param_dict = initialize_target_photometry_output_file(param_dict)
    for fits_filename in param_dict['fits_filenames']:
        param_dict = search_for_detection_in_extensions(param_dict,fits_filename)
        if param_dict['fits_filenames'][fits_filename]['target_ext_id'] != -1:
            if param_dict['fits_filenames'][fits_filename]['filter'] in filters:
                param_dict = perform_field_source_photometry(param_dict,fits_filename)
                param_dict = perform_target_source_photometry(param_dict,fits_filename)
                param_dict = write_target_photometry_to_file(param_dict,fits_filename)
            else:
                output_log_entry(param_dict['log_filepath'],'Skipping photometry of {:s} (non-standard filter).'.format(fits_filename))
        else:
            output_log_entry(param_dict['log_filepath'],'Skipping photometry of {:s} (object not in FOV).'.format(fits_filename))
        clean_output_files(param_dict,fits_filename)
    return param_dict

def convert_field_source_apert_params_to_pixels(param_dict):
    pixscale = param_dict['instrument_params']['pixscale']
    field_source_phot_radius_pix = param_dict['field_source_phot_radius_arcsec'] / pixscale
    sky_inner_r_pix              = param_dict['sky_inner_r_arcsec'] / pixscale
    sky_outer_r_pix              = param_dict['sky_outer_r_arcsec'] / pixscale
    param_dict['field_source_phot_radius_pix'] = field_source_phot_radius_pix
    param_dict['sky_inner_r_pix']              = sky_inner_r_pix
    param_dict['sky_outer_r_pix']              = sky_outer_r_pix    
    return param_dict

def convert_target_apert_params_to_pixels(param_dict):
    output_log_entry(param_dict['log_filepath'],'Converting target apertures from arcsec to pixels...')
    pixscale                 = param_dict['instrument_params']['pixscale']
    target_phot_radii_arcsec = param_dict['target_phot_radii_arcsec']
    target_phot_radii_pix    = [0 for idx in range(len(target_phot_radii_arcsec))]
    for idx in range(len(target_phot_radii_pix)):
        target_phot_radii_pix[idx] = target_phot_radii_arcsec[idx] / pixscale
    target_sky_inner_r_pix = param_dict['target_sky_inner_r_arcsec'] / pixscale
    target_sky_outer_r_pix = param_dict['target_sky_outer_r_arcsec'] / pixscale
    param_dict['target_phot_radii_pix']  = target_phot_radii_pix
    param_dict['target_sky_inner_r_pix'] = target_sky_inner_r_pix
    param_dict['target_sky_outer_r_pix'] = target_sky_outer_r_pix
    return param_dict

def search_for_detection_in_extensions(param_dict,fits_filename):
    pixscale      = param_dict['instrument_params']['pixscale']
    num_elements  = param_dict['instrument_params']['num_elements']
    first_element = param_dict['instrument_params']['first_element']
    output_log_entry(param_dict['log_filepath'],'Searching for target in {:s}...'.format(fits_filename))
    target_ext_id   = -1
    obj_ra  = param_dict['fits_filenames'][fits_filename]['ephem']['ra_deg']
    obj_dec = param_dict['fits_filenames'][fits_filename]['ephem']['dec_deg']
    for idx in range(0,num_elements):
        ext_id = first_element + idx
        ext_filename = fits_filename[:-5] + '_{:02d}.fits'.format(ext_id)
        if not os.path.exists(ext_filename):
            extract_extension_data(param_dict,fits_filename,ext_filename,ext_id)
            #with fits.open(fits_filename) as hdulist:
            #    hdr  = getheader(fits_filename,0)
            #    data = hdulist[0].data
            #    fits.writeto(ext_filename,data,hdr)
        ext_filename_wcs = compute_astrometric_solution(param_dict,ext_filename)
        if ext_filename_wcs != '':
            header = fits.getheader(ext_filename_wcs)
            w = WCS(header)
            with fits.open(ext_filename_wcs) as hdulist:
                data = hdulist[0].data
                npix_y,npix_x = data.shape
            obj_x,obj_y = w.wcs_world2pix(obj_ra,obj_dec,1)
            if obj_x >= 1 and obj_x <= npix_x and obj_y >= 1 and obj_y <= npix_y:
                output_log_entry(param_dict['log_filepath'],'Target in FOV of {:s} ({:.1f},{:.1f}).'.format(ext_filename_wcs,obj_x,obj_y))
                target_ext_id = ext_id
                break
            else:
                output_log_entry(param_dict['log_filepath'],'Target not in FOV of {:s} ({:.1f},{:.1f}).'.format(ext_filename_wcs,obj_x,obj_y))
                os.remove(ext_filename_wcs)
                obj_x,obj_y = 0,0
    if target_ext_id == -1:
        output_log_entry(param_dict['log_filepath'],'Target not found in {:s}.'.format(fits_filename))
    param_dict['fits_filenames'][fits_filename]['target_ext_id']    = target_ext_id
    param_dict['fits_filenames'][fits_filename]['ext_filename_wcs'] = ext_filename_wcs
    param_dict['fits_filenames'][fits_filename]['npix_x']           = npix_x
    param_dict['fits_filenames'][fits_filename]['npix_y']           = npix_y
    param_dict['fits_filenames'][fits_filename]['obj_x']            = obj_x
    param_dict['fits_filenames'][fits_filename]['obj_y']            = obj_y
    return param_dict


########## FIELD SOURCE PHOTOMETRY FUNCTIONS ##########

def perform_field_source_photometry(param_dict,fits_filename):
    output_log_entry(param_dict['log_filepath'],'Performing field source aperture photometry for {:s}...'.format(fits_filename))
    param_dict = find_field_sources(param_dict,fits_filename)
    param_dict = plot_field_source_apertures(param_dict,fits_filename)
    param_dict = measure_field_source_photometry(param_dict,fits_filename)
    param_dict = write_field_source_photometry_to_file(param_dict,fits_filename)
    param_dict = write_calibration_stars_to_file(param_dict,fits_filename)
    output_log_entry(param_dict['log_filepath'],'Field source aperture photometry for {:s} complete.'.format(fits_filename))
    return param_dict

def find_field_sources(param_dict,fits_filename):
    # Extracts field sources from image using DAOStarFinder
    # Based on code in https://photutils.readthedocs.io/en/stable/detection.html
    ext_filename_wcs       = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']
    ext_filename_wcs_bgsub = ext_filename_wcs[:-5] + '_bgsub.fits'
    param_dict['fits_filenames'][fits_filename]['ext_filename_wcs_bgsub'] = ext_filename_wcs_bgsub
    
    with fits.open(ext_filename_wcs) as hdu_data:
        data          = hdu_data[0].data
        hdr           = hdu_data[0].header
        sigma_clip    = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(data,(50,50),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
        data_bgsub    = data - bkg.background
        mean,median,std = sigma_clipped_stats(data_bgsub,sigma=3.0)
        fits.writeto(ext_filename_wcs_bgsub,data_bgsub,hdr,overwrite=True,checksum=True)
        
    #param_dict['fits_filenames'][fits_filename]['nx']     = nx
    #param_dict['fits_filenames'][fits_filename]['ny']     = ny
    #param_dict['fits_filenames'][fits_filename]['mean']   = mean
    #param_dict['fits_filenames'][fits_filename]['median'] = median
    #param_dict['fits_filenames'][fits_filename]['std']    = std
    
    daofind           = DAOStarFinder(fwhm=param_dict['daofind_fwhm'],threshold=param_dict['daofind_sigma_threshold']*std)
    field_sources     = daofind(data_bgsub)
    if field_sources is not None:
        num_field_sources = len(field_sources['xcentroid'])
    else:
        num_field_sources = 0
    
    if num_field_sources != 0:
        header = fits.getheader(ext_filename_wcs)
        w = WCS(header)
        field_sources['ra'],field_sources['dec'] = w.wcs_pix2world(field_sources['xcentroid'],field_sources['ycentroid'],1)
        field_sources['radec']                   = SkyCoord(field_sources['ra'],field_sources['dec'],unit='deg')
        field_sources['positions'] = np.transpose((field_sources['xcentroid'],field_sources['ycentroid']))

        output_filename = ext_filename_wcs[:-5]+'.field_sources.txt'
        with open(output_filename,'w') as of:
            of.write('star_id     RA (deg)   Dec (deg)     xcoord    ycoord\n')
            for idx in range(len(field_sources['ra'])):
                of.write('  {:03d}    {:11.7f} {:11.7f}  {:9.3f} {:9.3f}\n'.format(idx+1,field_sources['ra'][idx],field_sources['dec'][idx],field_sources['xcentroid'][idx],field_sources['ycentroid'][idx]))
    
    #print(field_sources.colnames)
    output_log_entry(param_dict['log_filepath'],'{:d} sources found'.format(num_field_sources))
    
    param_dict['fits_filenames'][fits_filename]['field_sources'] = field_sources
    
    return param_dict

#############################################################
## Define functions to plot source positions and apertures ##
#############################################################

def plot_field_source_apertures(param_dict,fits_filename):
    # Plots all source positions and apertures identified by DAOStarFinder
    output_log_entry(param_dict['log_filepath'],'Plotting field source photometry apertures for {:s}...'.format(fits_filename))
    ext_filename_wcs = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']

    data,error,imwcs = extract_image_file_data(ext_filename_wcs)
    ny,nx = data.shape
    
    field_source_phot_radius_pix = param_dict['field_source_phot_radius_pix']
    sky_inner_r_pix              = param_dict['sky_inner_r_pix']
    sky_outer_r_pix              = param_dict['sky_outer_r_pix']
    
    field_sources = param_dict['fits_filenames'][fits_filename]['field_sources']
    positions     = np.transpose((field_sources['xcentroid'],field_sources['ycentroid']))

    aperture         = CircularAperture(positions,r=field_source_phot_radius_pix)
    annulus_aperture = CircularAnnulus(positions,r_in=sky_inner_r_pix,r_out=sky_outer_r_pix)

    norm = simple_norm(data,'sqrt',percent=99)
    plt.imshow(data,norm=norm,interpolation='nearest')
    plt.xlim(0,nx)
    plt.ylim(0,ny)

    ap_patches  = aperture.plot(color='white',lw=1,label='Photometry aperture')
    ann_patches = annulus_aperture.plot(color='red',lw=1,label='Background annulus')
    handles     = (ap_patches[0],ann_patches[0])
    plt.legend(loc=(0.17,0.05),facecolor='#458989',labelcolor='white',handles=handles,prop={'weight':'bold','size':11})

    # Save to a File
    plot_filepath = ext_filename_wcs[:-5]+'.field_source_apertures.pdf'
    plt.savefig(plot_filepath,format = 'pdf', transparent=True)
    plt.clf()
    plt.cla()
    plt.close()

    output_log_entry(param_dict['log_filepath'],'Field source photometry apertures plotted ({:s}).'.format(plot_filepath))
    
    return param_dict


#####################################################
## Define functions to perform aperture photometry ##
#####################################################

def measure_field_source_photometry(param_dict,fits_filename):
    # photometry stored as Astropy tables in photometry dictionary indexed by image filename
    # aperture positions stored in param_dict['fits_filenames'][fits_filename]['field_sources']['positions']
    # photometry table columns: aperture_sum,aperture_sum_err,aperture_area,annulus_median,aperture_bkg,aperture_sum_bkgsub
    # "print(field_source_phot_table.colnames)" to access column names of astropy table

    field_source_phot_table = Table()
    ext_filename_wcs = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']

    # Set photometry aperture positions and sizes
    positions                       = param_dict['fits_filenames'][fits_filename]['field_sources']['positions']
    field_source_phot_radius_arcsec = param_dict['field_source_phot_radius_arcsec']
    field_source_phot_radius_pix    = param_dict['field_source_phot_radius_pix']
    aperture                        = CircularAperture(positions,r=field_source_phot_radius_pix)
    output_log_entry(param_dict['log_filepath'],'Photometry radius r = {:.2f} arcsec ({:.6f} px)'.format(field_source_phot_radius_arcsec,field_source_phot_radius_pix))
    
    # Read in image data
    data,error,imwcs = extract_image_file_data(ext_filename_wcs)
    
    # Compute image background properties from sky annuli
    sky_inner_r_pix  = param_dict['sky_inner_r_pix']
    sky_outer_r_pix  = param_dict['sky_outer_r_pix']
    annulus_aperture = CircularAnnulus(positions,r_in=sky_inner_r_pix,r_out=sky_outer_r_pix)
    annulus_mask     = annulus_aperture.to_mask(method='center')
    bkg_median,bkg_stdev = [],[]
    for mask in annulus_mask:
        annulus_data    = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _,median_sigclip,stdev_sigclip = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
        bkg_stdev.append(stdev_sigclip)
    bkg_median = np.array(bkg_median)
    bkg_stdev  = np.array(bkg_stdev)

    # Run aperture photometry function
    phot = aperture_photometry(data,aperture,method='exact',error=error)
    phot['aperture_bkg']        = bkg_median*aperture.area
    phot['aperture_sum_bkgsub'] = phot['aperture_sum'] - phot['aperture_bkg']

    # Compute zeropoint from field sources if image uses a standard broadband filter
    if param_dict['fits_filenames'][fits_filename]['filter'] in ['g_sdss','r_sdss','i_sdss','z_sdss','B','V','R','I']:
        param_dict = create_field_source_table(param_dict,fits_filename,phot)
        param_dict = compute_zero_point(param_dict,fits_filename)
        zeropoint  = ufloat(param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag'],param_dict['fits_filenames'][fits_filename]['median_zeropoint_err'])

    # Compute magnitudes of field sources from measured fluxes
    exptime = param_dict['fits_filenames'][fits_filename]['exposure_time']
    phot['mag'],phot['magerr'] = fluxes2mags(phot['aperture_sum_bkgsub'],phot['aperture_sum_err'],zeropoint,exptime)
    #phot['target_mag']    = [0.0 for idx in range(0,len(phot['mag']))]
    #phot['target_magerr'] = [0.0 for idx in range(0,len(phot['mag']))]
    #phot['mag_corr']      = [0.0 for idx in range(0,len(phot['mag']))]
    #phot['magerr_corr']   = [0.0 for idx in range(0,len(phot['mag']))]
    
    # Transfer data from aperture photometry output table to phot_table and convert to mag
    field_source_phot_table.add_column(phot['aperture_sum'],name='aperture_sum')
    field_source_phot_table.add_column(phot['aperture_sum_err'],name='aperture_sum_err')
    field_source_phot_table.add_column(aperture.area,name='aperture_area')
    field_source_phot_table.add_column(bkg_median,name='annulus_median')
    field_source_phot_table.add_column(bkg_stdev,name='annulus_stdev')
    field_source_phot_table.add_column(phot['aperture_bkg'],name='aperture_bkg')
    field_source_phot_table.add_column((phot['aperture_sum_bkgsub']),name='aperture_sum_bkgsub')
    field_source_phot_table.add_column((phot['mag']),name='mag')
    field_source_phot_table.add_column((phot['magerr']),name='magerr')
    #field_source_phot_table.add_column((phot['target_mag']),name='target_mag')
    #field_source_phot_table.add_column((phot['target_magerr']),name='target_magerr')
    #field_source_phot_table.add_column((phot['mag_corr']),name='mag_corr')
    #field_source_phot_table.add_column((phot['magerr_corr']),name='magerr_corr')
    
    param_dict['fits_filenames'][fits_filename]['field_source_phot_table'] = field_source_phot_table    

    return param_dict


def create_field_source_table(param_dict,fits_filename,phot):
    # Create field source table with aperture photometry
    field_source_table = Table()
    output_log_entry(param_dict['log_filepath'],'Creating field source table from param_dict...')    
    x     = param_dict['fits_filenames'][fits_filename]['field_sources']['xcentroid']
    y     = param_dict['fits_filenames'][fits_filename]['field_sources']['ycentroid']
    ra    = param_dict['fits_filenames'][fits_filename]['field_sources']['ra']
    dec   = param_dict['fits_filenames'][fits_filename]['field_sources']['dec']
    flux  = phot['aperture_sum_bkgsub']
    dflux = phot['aperture_sum_err']
    field_source_table = Table([x,y,ra,dec,flux,dflux], names=('x','y','ra','dec','flux','dflux'), dtype=('f8','f8','f8','f8','f8','f8'))
    #num_field_sources  = len(field_source_table)
    idx,field_sources_removed = 0,0
    while idx < len(field_source_table):
        # Remove stars with negative or zero fluxes or dflux > flux
        if field_source_table[idx]['flux'] <= 0 or field_source_table[idx]['flux'] < field_source_table[idx]['dflux']:
            field_source_table.remove_row(idx)
            field_sources_removed += 1
        else:
            idx += 1
    output_log_entry(param_dict['log_filepath'],'Removed {:d} bad sources (negative flux or large fluxerr) from field source table.'.format(field_sources_removed))
    field_source_table['SNR'] = -999.0
    num_sources = len(field_source_table)
    for idx in range(0,num_sources):
        field_source_table[idx]['SNR'] = field_source_table[idx]['flux'] / field_source_table[idx]['dflux']
    field_source_table['refcat_ra']   = -999.0
    field_source_table['refcat_dec']  = -999.0
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
    field_source_table['R_mag']       = -999.0
    field_source_table['R_err']       = -999.0
    field_source_table['I_mag']       = -999.0
    field_source_table['I_err']       = -999.0
    field_source_table['target_dist'] = -999.0
    field_source_table['zero_point']  = -999.0
    field_source_table['zpoint_err']  = -999.0
    param_dict['fits_filenames'][fits_filename]['field_source_table'] = field_source_table
    output_log_entry(param_dict['log_filepath'],'Field source table created.')
    return param_dict


########## TARGET PHOTOMETRY FUNCTIONS ##########

def perform_target_source_photometry(param_dict,fits_filename):
    param_dict = find_target_candidate_sources(param_dict,fits_filename)
    param_dict = find_target_source(param_dict,fits_filename)
    param_dict = plot_target_apertures(param_dict,fits_filename)
    param_dict = measure_target_photometry(param_dict,fits_filename)
    return param_dict

def initialize_target_photometry_output_file(param_dict):
    target_phot_radii_arcsec = param_dict['target_phot_radii_arcsec']
    for target_phot_radius_arcsec in target_phot_radii_arcsec:
        output_filename = 'target_phot_results.ap{:04.1f}arcsec.txt'.format(target_phot_radius_arcsec)
        target = param_dict['target_name']
        telescope = param_dict['telinst']
        with open(output_filename,'w') as of:
            of.write('target                telescope             filename                                                     jdmid  ')
            of.write('UTdate      UTtimemid     UTtimemidHr  exptime  ')
            of.write('filter  airmass      EphemRA     EphemDec  ')
            of.write('     RA_rate      Dec_rate       tot_rate  expt_trail_pa  expt_trail_len  heliodist   geodist  ')
            of.write('phsang    PA_AS   PA_NHV  orbplang  ')
            of.write('trueanom  Vmag  Nmag  Tmag  RA_3sigma  Dec_3sigma  ')
            of.write('   xcoord     ycoord             flux            dflux   ')
            of.write('     SNR           RA          Dec  zeropt  zpterr  ')
            of.write('nrefstars  apert     mag  magerr  limit_ps  limit_sb')
            of.write('\n')    
    return param_dict



##### FILE COMPRESSION FUNCTIONS #####

def fpack(filepath,config):
    filepath_fz = ''
    if 'cmd_fpack' in config:
        try:
            cmd = [config['cmd_fpack'],filepath]
            process = subprocess.call(cmd)
            os.remove(filepath)
            filepath_fz = filepath + '.fz'
        except:
            print('Function failed: fpack() (check if file is a FITS file)')
    else:
        print('Function failed: fpack() (local command not found)')
    return filepath_fz

def funpack(filepath_fz,config):
    filepath = ''
    if 'cmd_funpack' in config:
        try:
            cmd = [config['cmd_funpack'],filepath_fz]
            process = subprocess.call(cmd)
            os.remove(filepath_fz)
            filepath = filepath_fz[:-3]
        except:
            print('Function failed: funpack() (check if file is a fpacked FITS file)')
    else:
        print('Function failed: funpack() (local command not found)')
    return filepath


##### UTILITY FUNCTIONS #####

def file_len(fname):
    p = subprocess.Popen(['wc','-l',fname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    result,err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def remove_extra_characters(string):
    new_string = string.replace('-','')
    new_string = new_string.replace('/','')
    new_string = new_string.replace('+','')
    new_string = new_string.replace(' ','')
    return new_string


##### SETUP FUNCTIONS #####


##### LOGGING FUNCTIONS #####


def output_log_entry(log_filepath,log_entry):
    with open(log_filepath,'a') as log_file:
        print('{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
        log_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
    return None






##### SSOIS SEARCH FUNCTIONS #####





##### ASTROMETRY FUNCTIONS #####

def compute_astrometric_solution(param_dict,fits_filename):
    wcs_filename = fits_filename[:-5] + '_wcs.fits'
    if not os.path.exists(wcs_filename):
        output_log_entry(param_dict['log_filepath'],'Computing astrometric solution for {:s}...'.format(fits_filename))
        pixscale = param_dict['instrument_params']['pixscale']    
        if pixscale != 0:
            cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7 = 'solve-field','--scale-units','arcsecperpix','--scale-low','{:.3f}'.format(pixscale-0.1),'--scale-high','{:.3f}'.format(pixscale+0.1)
            cmd = [cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,fits_filename]
            process = subprocess.call(cmd)        
        else:
            cmd1,cmd2,cmd3 = 'solve-field','--scale-low','1'
            cmd = [cmd1,cmd2,cmd3,fits_filename]
            process = subprocess.call(cmd)
        if os.path.isfile(fits_filename[:-5] + '.solved'):
            wcs_filename = fits_filename[:-5] + '_wcs.fits'
            os.rename(fits_filename[:-5]+'.new',wcs_filename)
            output_log_entry(param_dict['log_filepath'],'Astrometric solution successfully completed.')
            clean_wcs_output(fits_filename,param_dict)
            os.remove(fits_filename)
        else:
            wcs_filename = ''
            output_log_entry(param_dict['log_filepath'],'Astrometric solution failed.')
    else:
        output_log_entry(param_dict['log_filepath'],'Astrometric solution for {:s} already computed.'.format(fits_filename))
        
    return wcs_filename

def clean_wcs_output(fits_filename,param_dict):
    output_log_entry(param_dict['log_filepath'],'Cleaning up astrometry.net output...')
    if os.path.exists(fits_filename[:-5]+'-indx.png'):  os.remove(fits_filename[:-5]+'-indx.png')
    if os.path.exists(fits_filename[:-5]+'-ngc.png'):   os.remove(fits_filename[:-5]+'-ngc.png')
    if os.path.exists(fits_filename[:-5]+'-objs.png'):  os.remove(fits_filename[:-5]+'-objs.png')
    if os.path.exists(fits_filename[:-5]+'-indx.xyls'): os.remove(fits_filename[:-5]+'-indx.xyls')
    if os.path.exists(fits_filename[:-5]+'.axy'):       os.remove(fits_filename[:-5]+'.axy')
    if os.path.exists(fits_filename[:-5]+'.corr'):      os.remove(fits_filename[:-5]+'.corr')
    if os.path.exists(fits_filename[:-5]+'.match'):     os.remove(fits_filename[:-5]+'.match')
    if os.path.exists(fits_filename[:-5]+'.rdls'):      os.remove(fits_filename[:-5]+'.rdls')
    if os.path.exists(fits_filename[:-5]+'.solved'):    os.remove(fits_filename[:-5]+'.solved')
    if os.path.exists(fits_filename[:-5]+'.wcs'):       os.remove(fits_filename[:-5]+'.wcs')
    return None



def extract_image_file_data(fits_filename):
    with fits.open(fits_filename) as hdu:
        data_tmp = hdu[0].data
        data     = ma.masked_invalid(data_tmp) # mask NaN values
        ny,nx    = data.shape
        imwcs    = None
        mean,median,std = sigma_clipped_stats(data,sigma=3.0)
        error    = [[std for idx in range(nx)] for idx in range(ny)]
    return data,error,imwcs

def extract_extension_data(param_dict,fits_filename,ext_filename,ext_id):
    if param_dict['telinst'] == 'SDSS':
        with fits.open(fits_filename) as hdulist:
            hdr  = getheader(fits_filename,ext_id)
            data = hdulist[ext_id].data
            fits.writeto(ext_filename,data,hdr)
    if param_dict['telinst'] == 'CTIO-4m/DECam':
        with fits.open(fits_filename) as hdulist:
            hdr0  = getheader(fits_filename,0)
            hdr1  = getheader(fits_filename,ext_id)
            data = hdulist[ext_id].data
            fits.writeto(ext_filename,data,hdr0+hdr1)
    if param_dict['telinst'] == 'CFHT/MegaCam':
        with fits.open(fits_filename) as hdulist:
            hdr0  = getheader(fits_filename,0)
            hdr1  = getheader(fits_filename,ext_id)
            data = hdulist[ext_id].data
            fits.writeto(ext_filename,data,hdr0+hdr1)
    return None


########## EPHEMERIS RETRIEVAL FUNCTIONS ##########

def run_horizons_queries_prelim(param_dict):
    target_name = param_dict['target_name']
    obs_code    = param_dict['instrument_params']['obs_code']  # observatory code

    query_results_ephems_filepath = param_dict['query_results_filepath'][:-4] + '_ephems.txt'

    with open(query_results_ephems_filepath,'w') as of:
        of.write('Object        Date       Time          JD              Tel/Inst               Exptime  Filter          Rdist     Ddist  PhsA  TrAnm  Vmag  Tmag  Nmag  Target                          PsAng  PsAMV  OrbPl  RA        Dec        RA          Dec         RA_rate  Dec_rate  RA_3sig  Dec_3sig  EclLon  EclLat  GlxLon  GlxLat\n')

    with open(param_dict['query_results_filepath'],'r') as query_results_file:
        for _ in range(1): #skip first header line
            next(query_results_file)
        for line in query_results_file:
            ephem       = {}
            line_data   = re.split(r'\t+',line.rstrip('\t'))
            image_name  = line_data[0].split()[0].strip()
            mjd         = float(str(line_data[1]))
            obs_time    = Time('{:f}'.format(mjd),format='mjd')
            date        = obs_time.iso[:10]
            time        = obs_time.iso[11:22]            
            filter_name = line_data[2]
            exptime     = float(str(line_data[3]))
            obj_ra_deg  = float(str(line_data[4]))
            obj_dec_deg = float(str(line_data[5]))
            obs_target  = line_data[6]

            obj = Horizons(id=param_dict['target_name'],location=obs_code,epochs=obs_time.jd)
            query_successful = True
            max_retries = 50
            output_log_entry(param_dict['log_filepath'],'Trying Horizons in comet mode...')
            for _ in range(max_retries):
                try:
                    if (not os.environ.get('PYTHONHTTPSVERIFY','') and getattr(ssl,'_create_unverified_context',None)):
                        ssl._create_default_https_context = ssl._create_unverified_context
                    ephems = obj.ephemerides(closest_apparition=True,no_fragments=True)
                    print(obj.uri)
                except Exception:
                    print('.',end='')
                else:
                    print('')
                    break
            else:
                print('\n')
                output_log_entry(param_dict['log_filepath'],'Horizons query (comet mode) failed. Maximum retries reached.')
                query_successful = False
        
            # try asteroid mode
            if not query_successful:
                query_successful = True
                output_log_entry(param_dict['log_filepath'],'Trying Horizons in asteroid mode...')
                for _ in range(max_retries):
                    try:
                        if (not os.environ.get('PYTHONHTTPSVERIFY','') and getattr(ssl,'_create_unverified_context',None)):
                            ssl._create_default_https_context = ssl._create_unverified_context
                        ephems = obj.ephemerides()
                    except Exception:
                        print('.',end='')
                    else:
                        print('')
                        break
                else:
                    print('\n')
                    output_log_entry(param_dict['log_filepath'],'Horizons query (asteroid mode) failed. Maximum retries reached.')
                    query_successful = False
            
            if query_successful:   
                #print(ephems.columns)
                ephem['v_mag'],ephem['t_mag'],ephem['n_mag']  = 99.9,99.9,99.9
                ephem['ra_deg']             = float(str(ephems[0]['RA']))
                ephem['dec_deg']            = float(str(ephems[0]['DEC']))
                ra_angle                    = Angle(ephem['ra_deg'] * u.deg)
                dec_angle                   = Angle(ephem['dec_deg'] * u.deg)
                ra_hms                      = ra_angle.to_string(unit=u.hour,sep=':',precision=0,pad=True)
                dec_dms                     = dec_angle.to_string(unit=u.degree,sep=':',precision=0,pad=True,alwayssign=True)
                ephem['ra_rate']            = float(str(ephems[0]['RA_rate']))    # arcsec/hr
                ephem['dec_rate']           = float(str(ephems[0]['DEC_rate']))   # arcsec/hr
                ephem['phase_angle']        = float(str(ephems[0]['alpha']))
                ephem['rdist']              = float(str(ephems[0]['r']))
                ephem['ddist']              = float(str(ephems[0]['delta']))
                if 'V' in ephems.columns:
                    try:
                        ephem['v_mag'] = float(str(ephems[0]['V']))
                    except ValueError:
                        ephem['v_mag'] = 99.9
                if 'Tmag' in ephems.columns:
                    try:
                        ephem['t_mag'] = float(str(ephems[0]['Tmag']))
                    except ValueError:
                        ephem['t_mag'] = 99.9
                if 'Nmag' in ephems.columns:
                    try:
                        ephem['n_mag'] = float(str(ephems[0]['Nmag']))
                    except ValueError:
                        ephem['n_mag'] = 99.9
                if ephem['v_mag'] == 99.9:
                    ephem['v_mag_str'] = ' -- '
                else:
                    ephem['v_mag_str'] = '{:4.1f}'.format(ephem['v_mag'])
                if ephem['t_mag'] == 99.9:
                    ephem['t_mag_str'] = ' -- '
                else:
                    ephem['t_mag_str'] = '{:4.1f}'.format(ephem['t_mag'])
                if ephem['n_mag'] == 99.9:
                    ephem['n_mag_str'] = ' -- '
                else:
                    ephem['n_mag_str'] = '{:4.1f}'.format(ephem['n_mag'])
                ephem['pa_as']               = float(str(ephems[0]['sunTargetPA']))
                ephem['pa_nhv']              = float(str(ephems[0]['velocityPA']))
                ephem['ecliptic_longitude']  = float(str(ephems[0]['EclLon']))
                ephem['ecliptic_latitude']   = float(str(ephems[0]['EclLat']))
                ephem['true_anomaly']        = float(str(ephems[0]['true_anom']))
                ephem['galactic_longitude']  = float(str(ephems[0]['GlxLon']))
                ephem['galactic_latitude']   = float(str(ephems[0]['GlxLat']))
                ephem['ra_3sigma']           = float(str(ephems[0]['RA_3sigma']))
                ephem['dec_3sigma']          = float(str(ephems[0]['DEC_3sigma']))
                ephem['orb_pl_angle']        = float(str(ephems[0]['OrbPlaneAng']))
                
                with open(query_results_ephems_filepath,'a') as of:
                    #of.write('Object        Date       Time          JD              Tel/Inst               Exptime  Filter          Rdist     Ddist  PhsA  TrAnm  Vmag  Tmag  Nmag    Target                          PsAng  PsAMV  OrbPl  RA        Dec        RA          Dec         RA_rate  Dec_rate  RA3sigma  DEC3sigma  EclLon  EclLat  GlxLon  GlxLat\n')
                    of.write('{:>12s}  {:<10s} {:<12s}  {:14.6f}  {:<21s}  {:7.1f}  {:<11s}  {:8.3f}  {:8.3f}  {:4.1f}  {:5.1f}  {:s}  {:s}  {:s}  {:<30s}  {:5.1f}  {:5.1f}  {:5.1f}  {:s}  {:s}  {:10.6f}  {:10.6f}  {:7.3f}  {:8.3f}  {:7.1f}  {:8.1f}  {:6.1f}  {:6.1f}  {:6.1f}  {:6.1f}\n'.format(target_name,date,time,obs_time.jd,param_dict['telinst'],exptime,filter_name,ephem['rdist'],ephem['ddist'],ephem['phase_angle'],ephem['true_anomaly'],ephem['v_mag_str'],ephem['t_mag_str'],ephem['n_mag_str'],obs_target,ephem['pa_as'],ephem['pa_nhv'],ephem['orb_pl_angle'],ra_hms,dec_dms,ephem['ra_deg'],ephem['dec_deg'],ephem['ra_rate'],ephem['dec_rate'],ephem['ra_3sigma'],ephem['dec_3sigma'],ephem['ecliptic_longitude'],ephem['ecliptic_latitude'],ephem['galactic_longitude'],ephem['galactic_latitude']))

    return param_dict


def run_horizons_queries(param_dict):
    target_name = param_dict['target_name']

    for fits_filename in param_dict['fits_filenames']:
        ephem = {}
        output_log_entry(param_dict['log_filepath'],'Running Horizons query for {:s} for {:s}...'.format(target_name,fits_filename))
        ra_deg,dec_deg,ra_rate,dec_rate,v_mag_str,t_mag_str,n_mag_str,ra_3sigma,dec_3sigma = 0,0,0,0,0,0,0,0,0

        obs_code    = param_dict['instrument_params']['obs_code']  # observatory code
        exptime     = param_dict['fits_filenames'][fits_filename]['exposure_time']
        filter_name = param_dict['fits_filenames'][fits_filename]['filter']
        jd_mid      = param_dict['fits_filenames'][fits_filename]['time_mid_jd']
        
        obj = Horizons(id=param_dict['target_name'],location=obs_code,epochs=jd_mid)
        query_successful = True
        max_retries = 50
        output_log_entry(param_dict['log_filepath'],'Trying Horizons in comet mode...')
        for _ in range(max_retries):
            try:
                if (not os.environ.get('PYTHONHTTPSVERIFY','') and getattr(ssl,'_create_unverified_context',None)):
                    ssl._create_default_https_context = ssl._create_unverified_context
                ephems = obj.ephemerides(closest_apparition=True,no_fragments=True)
                print(obj.uri)
            except Exception:
                print('.',end='')
            else:
                print('')
                break
        else:
            print('\n')
            output_log_entry(param_dict['log_filepath'],'Horizons query (comet mode) failed. Maximum retries reached.')
            query_successful = False
        
        # try asteroid mode
        if not query_successful:
            query_successful = True
            output_log_entry(param_dict['log_filepath'],'Trying Horizons in asteroid mode...')
            for _ in range(max_retries):
                try:
                    if (not os.environ.get('PYTHONHTTPSVERIFY','') and getattr(ssl,'_create_unverified_context',None)):
                        ssl._create_default_https_context = ssl._create_unverified_context
                    ephems = obj.ephemerides()
                except Exception:
                    print('.',end='')
                else:
                    print('')
                    break
            else:
                print('\n')
                output_log_entry(param_dict['log_filepath'],'Horizons query (asteroid mode) failed. Maximum retries reached.')
                query_successful = False
            
        if query_successful:   
            #print(ephems.columns)
            ephem['v_mag'],ephem['t_mag'],ephem['n_mag']  = 99.9,99.9,99.9
            ephem['ra_deg']             = float(str(ephems[0]['RA']))
            ephem['dec_deg']            = float(str(ephems[0]['DEC']))
            ra_angle                    = Angle(ephem['ra_deg'] * u.deg)
            dec_angle                   = Angle(ephem['dec_deg'] * u.deg)
            ra_hms                      = ra_angle.to_string(unit=u.hour,sep=':',precision=0,pad=True)
            dec_dms                     = dec_angle.to_string(unit=u.degree,sep=':',precision=0,pad=True,alwayssign=True)
            ephem['ra_rate']            = float(str(ephems[0]['RA_rate']))    # arcsec/hr
            ephem['dec_rate']           = float(str(ephems[0]['DEC_rate']))   # arcsec/hr
            ephem['airmass']            = float(str(ephems[0]['airmass']))
            ephem['tot_rate']           = ((ephem['ra_rate']**2)+(ephem['dec_rate']**2))**0.5
            ephem['trail_pa_expected']  = np.arctan(ephem['ra_rate']/ephem['dec_rate'])/math.pi*180
            ephem['trail_len_expected'] = ephem['tot_rate'] / 3600 * param_dict['fits_filenames'][fits_filename]['exposure_time']
            ephem['phase_angle']        = float(str(ephems[0]['alpha']))
            ephem['rdist']              = float(str(ephems[0]['r']))
            ephem['ddist']              = float(str(ephems[0]['delta']))
            if 'V' in ephems.columns:
                try:
                    ephem['v_mag'] = float(str(ephems[0]['V']))
                except ValueError:
                    ephem['v_mag'] = 99.9
                #if ephems[0]['V'].isnumeric(): v_mag = float(ephems[0]['V'])
            if 'Tmag' in ephems.columns:
                try:
                    ephem['t_mag'] = float(str(ephems[0]['Tmag']))
                except ValueError:
                    ephem['t_mag'] = 99.9
            if 'Nmag' in ephems.columns:
                try:
                    ephem['n_mag'] = float(str(ephems[0]['Nmag']))
                except ValueError:
                    ephem['n_mag'] = 99.9
            if ephem['v_mag'] == 99.9:
                ephem['v_mag_str'] = ' -- '
            else:
                ephem['v_mag_str'] = '{:4.1f}'.format(ephem['v_mag'])
            if ephem['t_mag'] == 99.9:
                ephem['t_mag_str'] = ' -- '
            else:
                ephem['t_mag_str'] = '{:4.1f}'.format(ephem['t_mag'])
            if ephem['n_mag'] == 99.9:
                ephem['n_mag_str'] = ' -- '
            else:
                ephem['n_mag_str'] = '{:4.1f}'.format(ephem['n_mag'])
            ephem['pa_as']               = float(str(ephems[0]['sunTargetPA']))
            ephem['pa_nhv']              = float(str(ephems[0]['velocityPA']))
            ephem['ecliptic_longitude']  = float(str(ephems[0]['EclLon']))
            ephem['ecliptic_latitude']   = float(str(ephems[0]['EclLat']))
            ephem['true_anomaly']        = float(str(ephems[0]['true_anom']))
            ephem['galactic_longitude']  = float(str(ephems[0]['GlxLon']))
            ephem['galactic_latitude']   = float(str(ephems[0]['GlxLat']))
            ephem['ra_3sigma']           = float(str(ephems[0]['RA_3sigma']))
            ephem['dec_3sigma']          = float(str(ephems[0]['DEC_3sigma']))
            ephem['orb_pl_angle']        = float(str(ephems[0]['OrbPlaneAng']))
            #ephem['solar_elongation']   = float(str(ephems[0]['elong']))
            #lunar_illumination          = float(str(ephems[0]['lunar_illum']))
            #lunar_elongation            = float(str(ephems[0]['lunar_elong']))


        param_dict['fits_filenames'][fits_filename]['ephem'] = ephem
    
    return param_dict






########## APERTURE SIZE CONVERSION FUNCTIONS ##########





##### PS1 MAGNITUDE TRANSFORMATION FUNCTION #####

def convert_ps1_mags(gp1_mag,gp1_err,rp1_mag,rp1_err,ip1_mag,ip1_err,zp1_mag,zp1_err):
    # Convert PS1 magnitudes to SDSS and Johnson-Cousins systems
    magnitudes = {'g_sdss':{'mag':0,'err':0},'r_sdss':{'mag':0,'err':0},
                  'i_sdss':{'mag':0,'err':0},'z_sdss':{'mag':0,'err':0},
                  'B':{'mag':0,'err':0},'V':{'mag':0,'err':0},'R':{'mag':0,'err':0},'I':{'mag':0,'err':0}}

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
    
    magnitudes['g_sdss']['mag'],magnitudes['g_sdss']['err'] = magavg([g_sdss_mag1,g_sdss_mag2],[gp1_err,gp1_err])
    magnitudes['r_sdss']['mag'],magnitudes['r_sdss']['err'] = magavg([r_sdss_mag1,r_sdss_mag2],[rp1_err,rp1_err])
    magnitudes['i_sdss']['mag'],magnitudes['i_sdss']['err'] = magavg([i_sdss_mag1,i_sdss_mag2],[ip1_err,ip1_err])
    magnitudes['z_sdss']['mag'],magnitudes['z_sdss']['err'] = magavg([z_sdss_mag1,z_sdss_mag2],[zp1_err,zp1_err])
    magnitudes['B']['mag'],     magnitudes['B']['err']      = magavg([B_mag1, B_mag2] ,[gp1_err,gp1_err])
    magnitudes['V']['mag'],     magnitudes['V']['err']      = magavg([V_mag1, V_mag2] ,[rp1_err,rp1_err])
    magnitudes['R']['mag'],     magnitudes['R']['err']      = magavg([Rc_mag1,Rc_mag2],[rp1_err,rp1_err])
    magnitudes['I']['mag'],     magnitudes['I']['err']      = magavg([Ic_mag1,Ic_mag2],[ip1_err,ip1_err])
    
    return magnitudes


##### MEAN/MEDIAN MAGNITUDE CALCULATION FUNCTIONS #####

def magavg(mag_array,magerr_array):
    avgmag,avgmagerr = -999,-999
    num_mags = len(mag_array)
    intensity_total = 0
    error_total = 0
    # initialize intensity array
    intensity = [[0 for idx1 in range(2)] for idx2 in range(num_mags)]
    for idx in range(0,num_mags):
        intensity[idx][0] = 10**(0.4*(0.-mag_array[idx]))
        intensity[idx][1] = ((magerr_array[idx]**2) * ((-0.921*(10**(-0.4*mag_array[idx])))**2.))**0.5
        intensity_total   = intensity_total + intensity[idx][0] / (intensity[idx][1]**2)
        error_total       = error_total + (intensity[idx][1]**2)**(-1)
    avgintensity = intensity_total / error_total
    avgintensityerr = error_total ** (-0.5)
    avgmag = -2.5 * math.log10(avgintensity)
    avgmagerr = 1.08574 * avgintensityerr/avgintensity
    return avgmag,avgmagerr


def magmedian(mag_array):
    flux_array      = [0 for idx in range(len(mag_array))]
    flux_dev2_array = [0 for idx in range(len(mag_array))]
    mag_dev2_array  = [0 for idx in range(len(mag_array))]
    for idx in range(len(mag_array)):
        flux_array[idx] = 10**(0.4*(0.-mag_array[idx]))
    median_flux = statistics.median(flux_array)
    median_mag = -2.5 * math.log10(median_flux)
    for idx in range(len(flux_array)):
        flux_dev2_array[idx] = (flux_array[idx] - median_flux)**2
        #mag_dev2_array[idx]  = (mag_array[idx] - median_mag)**2
    flux_stddev = (np.sum(flux_dev2_array)/len(flux_dev2_array))**0.5
    #mag_stddev  = (np.sum(mag_dev2_array)/len(mag_dev2_array))**0.5
    mag_stddev_fluxspace_dn = median_mag - (-2.5 * math.log10(median_flux+flux_stddev))
    if (median_flux-flux_stddev) > 0:
        mag_stddev_fluxspace_up = (-2.5 * math.log10(median_flux-flux_stddev)) - median_mag
    else:
        mag_stddev_fluxspace_up = median_mag
    median_mag_err = max([mag_stddev_fluxspace_dn,mag_stddev_fluxspace_up])
    return median_mag,median_mag_err


def fluxes2mags(flux,fluxerr,zeropoint,exptime):
    mag    = [0 for idx in range(len(flux))]
    magerr = [0 for idx in range(len(flux))]
    for idx in range(len(flux)):
        if flux[idx] > 0 and fluxerr[idx] > 0:
            if isinstance(zeropoint,uncertainties.UFloat):
                mag_ufloat  = zeropoint.n - 2.5 * unumpy.log10(ufloat(flux[idx],fluxerr[idx])/exptime)
            else:
                mag_ufloat  = zeropoint - 2.5 * unumpy.log10(ufloat(flux[idx],fluxerr[idx])/exptime)
            mag[idx]    = mag_ufloat.n
            magerr[idx] = mag_ufloat.s
            #magerr[idx] = 1.08574 * fluxerr[idx]/flux[idx]
            #print('{:.3f}+/-{:.3f}'.format(mag[idx],magerr[idx]))
        else:
            mag[idx]    = 99.999
            magerr[idx] = 9.999
    return mag,magerr


##### FIELD SOURCE PHOTOMETRY FUNCTIONS #####



def compute_zero_point(param_dict,fits_filename):
    # Compute zero points
    field_source_table = param_dict['fits_filenames'][fits_filename]['field_source_table']
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs_bgsub']
    
    filter_name = param_dict['fits_filenames'][fits_filename]['filter']
    exptime     = param_dict['fits_filenames'][fits_filename]['exposure_time']
    x           = field_source_table['x']
    y           = field_source_table['y']
    ra          = field_source_table['ra']
    dec         = field_source_table['dec']
    flux        = field_source_table['flux']
    dflux       = field_source_table['dflux']
    SNR         = flux/dflux
    
    calib_star_table = Table(names=('x','y','ra','dec','flux','dflux','SNR','refcat_ra','refcat_dec',\
                                    'gp1_mag','gp1_err','rp1_mag','rp1_err',\
                                    'ip1_mag','ip1_err','zp1_mag','zp1_err',\
                                    'g_sdss_mag','g_sdss_err','r_sdss_mag','r_sdss_err',\
                                    'i_sdss_mag','i_sdss_err','z_sdss_mag','z_sdss_err',\
                                    'B_mag','B_err','V_mag','V_err','R_mag','R_err','I_mag','I_err',\
                                    'target_dist','zero_point','zpoint_err'))

    for idx in range(0,len(field_source_table)):
        refcat_photometry_dict = find_individual_refcat_source_photometry(ra[idx],dec[idx],param_dict)
        if refcat_photometry_dict is not None:
            #print(refcat_photometry_dict)
            targ_magerr  = ufloat(refcat_photometry_dict[filter_name]['mag'],refcat_photometry_dict[filter_name]['err'])
            star_fluxerr = ufloat(flux[idx],dflux[idx])
            zero_point   = targ_magerr + 2.5 * unumpy.log10(star_fluxerr/exptime)
    
            refcat_ra             = refcat_photometry_dict['refcat_ra']
            refcat_dec            = refcat_photometry_dict['refcat_dec']
            gp1_mag,gp1_err       = refcat_photometry_dict['gp1']['mag'],   refcat_photometry_dict['gp1']['err']
            rp1_mag,rp1_err       = refcat_photometry_dict['rp1']['mag'],   refcat_photometry_dict['rp1']['err']
            ip1_mag,ip1_err       = refcat_photometry_dict['ip1']['mag'],   refcat_photometry_dict['ip1']['err']
            zp1_mag,zp1_err       = refcat_photometry_dict['zp1']['mag'],   refcat_photometry_dict['zp1']['err']
            g_sdss_mag,g_sdss_err = refcat_photometry_dict['g_sdss']['mag'],refcat_photometry_dict['g_sdss']['err']
            r_sdss_mag,r_sdss_err = refcat_photometry_dict['r_sdss']['mag'],refcat_photometry_dict['r_sdss']['err']
            i_sdss_mag,i_sdss_err = refcat_photometry_dict['i_sdss']['mag'],refcat_photometry_dict['i_sdss']['err']
            z_sdss_mag,z_sdss_err = refcat_photometry_dict['z_sdss']['mag'],refcat_photometry_dict['z_sdss']['err']
            B_mag,B_err           = refcat_photometry_dict['B']['mag'],     refcat_photometry_dict['B']['err']
            V_mag,V_err           = refcat_photometry_dict['V']['mag'],     refcat_photometry_dict['V']['err']
            R_mag,R_err           = refcat_photometry_dict['R']['mag'],     refcat_photometry_dict['R']['err']
            I_mag,I_err           = refcat_photometry_dict['I']['mag'],     refcat_photometry_dict['I']['err']
            
            target_dist = refcat_photometry_dict['target_dist']
            
            calib_star_table.add_row((x[idx],y[idx],ra[idx],dec[idx],flux[idx],dflux[idx],SNR[idx],refcat_ra,refcat_dec,
                                      gp1_mag,gp1_err,rp1_mag,rp1_err,ip1_mag,ip1_err,zp1_mag,zp1_err,
                                      g_sdss_mag,g_sdss_err,r_sdss_mag,r_sdss_err,i_sdss_mag,i_sdss_err,z_sdss_mag,z_sdss_err,
                                      B_mag,B_err,V_mag,V_err,R_mag,R_err,I_mag,I_err,
                                      target_dist,zero_point.n,zero_point.s))
    
    output_log_entry(param_dict['log_filepath'],'Removing sources with outlier zero point magnitudes...')
    zeropoints = calib_star_table['zero_point']
    #print(zeropoints)
    mu,sigma = generate_zeropoint_histogram(calib_star_table['zero_point'],fits_filename,'histogram1',param_dict)
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
    output_log_entry(param_dict['log_filepath'],'{:d} outlier zeropoints removed, {:d} zeropoints remaining.'.format(zpoints_removed,(num_calib_stars-zpoints_removed)))
    
    if (num_calib_stars-zpoints_removed) < 3:
        output_log_entry(param_dict['log_filepath'],'Insufficient zeropoints (n<3) remaining for compute_zero_point() to continue for {:s}'.format(ext_filename_wcs))
        median_zeropoint_mag,median_zeropoint_err = 0,0
    else:
        mu,sigma = generate_zeropoint_histogram(calib_star_table['zero_point'],fits_filename,'histogram2',param_dict)
        median_zeropoint_mag,median_zeropoint_err = magmedian(calib_star_table['zero_point'])
        output_log_entry(param_dict['log_filepath'],'Number of stars used for photometric solution: {:d}'.format(num_calib_stars))
        output_log_entry(param_dict['log_filepath'],'Initial zero point median: {:.3f}'.format(zeropoint_median))
        output_log_entry(param_dict['log_filepath'],'Final median zero point: {:.3f}+/-{:.3f}'.format(median_zeropoint_mag,median_zeropoint_err))
    
    param_dict['fits_filenames'][fits_filename]['calib_star_table'] = calib_star_table
    param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag'] = median_zeropoint_mag
    param_dict['fits_filenames'][fits_filename]['median_zeropoint_err'] = median_zeropoint_err
    
    return param_dict


def find_individual_refcat_source_photometry(source_ra,source_dec,param_dict):
    refcat_photometry_dict      = None
    search_radius_arcsec        = 10
    cmd_refcat                  = param_dict['config']['cmd_refcat']
    refcat_sources_filename_tmp = 'refcat_sources_temp.txt'
    
    with open(refcat_sources_filename_tmp,'w') as outf:
        cmd_ra      = '{:f}'.format(source_ra)
        cmd_dec     = '{:f}'.format(source_dec)
        cmd_radius  = '{:f}'.format(search_radius_arcsec/3600)
        directories = '{:s},{:s},{:s},{:s},{:s}'.format(param_dict['config']['refcat_dir_00_m_16'],
                        param_dict['config']['refcat_dir_16_m_17'],param_dict['config']['refcat_dir_17_m_18'],
                        param_dict['config']['refcat_dir_18_m_19'],param_dict['config']['refcat_dir_19_m_20'])
        #print(directories)
        var_fields = 'RA,Dec,g,dg,r,dr,i,di,z,dz'
        cmd_dir,cmd_rad,cmd_var,cmd_nohdr,cmd_bin = '-dir','-rad','-var','-nohdr','-bin'
        cmd = [cmd_refcat,cmd_ra,cmd_dec,cmd_dir,directories,cmd_bin,cmd_rad,cmd_radius,cmd_var,var_fields,cmd_nohdr]
        #print(cmd)
        subprocess.call(cmd,stdout=outf)

    # parse local file into astropy.table object
    if file_len(refcat_sources_filename_tmp) > 0:
        objRA  = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(0))
        objDec = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(1))
        gmag   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(2))
        gerr   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(3))
        rmag   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(4))
        rerr   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(5))
        imag   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(6))
        ierr   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(7))
        zmag   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(8))
        zerr   = np.genfromtxt(refcat_sources_filename_tmp,skip_header=0,usecols=(9))

        if objRA.size == 1:
            objRA  = np.array([float(objRA)])
            objDec = np.array([float(objDec)])
            gmag   = np.array([float(gmag)])
            gerr   = np.array([float(gerr)])
            rmag   = np.array([float(rmag)])
            rerr   = np.array([float(rerr)])
            imag   = np.array([float(imag)])
            ierr   = np.array([float(ierr)])
            zmag   = np.array([float(zmag)])
            zerr   = np.array([float(zerr)])

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

        # Find matching refcat source
        idx,matched_refcat_ra,matched_refcat_dec = 0,0,0
        star_matched = False
        while not star_matched and idx < num_refcat_sources:
            refcat_ra  = refcat_table[idx]['right_ascension']
            refcat_dec = refcat_table[idx]['declination']
            c1         = SkyCoord(source_ra,source_dec,unit='deg',frame='icrs')
            c2         = SkyCoord(refcat_ra,refcat_dec,unit='deg',frame='icrs')
            sep        = c1.separation(c2)
            if sep.arcsecond < 2.6:  # If target matches within 2.5 arcsec
                star_matched = True
                matched_refcat_ra  = refcat_ra
                matched_refcat_dec = refcat_dec
                target_dist        = sep.arcsecond
                gp1_mag            = refcat_table[idx]['gp1_mag']
                gp1_err            = refcat_table[idx]['gp1_err']
                rp1_mag            = refcat_table[idx]['rp1_mag']
                rp1_err            = refcat_table[idx]['rp1_err']
                ip1_mag            = refcat_table[idx]['ip1_mag']
                ip1_err            = refcat_table[idx]['ip1_err']
                zp1_mag            = refcat_table[idx]['zp1_mag']
                zp1_err            = refcat_table[idx]['zp1_err']
                output_log_entry(param_dict['log_filepath'],'Field source at {:.6f},{:.6f} --> refcat source at {:.6f},{:.6f}'.format(source_ra,source_dec,refcat_ra,refcat_dec))
            idx += 1

        if star_matched:
            magnitudes = convert_ps1_mags(gp1_mag,gp1_err,rp1_mag,rp1_err,ip1_mag,ip1_err,zp1_mag,zp1_err)
            refcat_photometry_dict = {'source_ra':source_ra,'source_dec':source_dec,
                                      'refcat_ra':matched_refcat_ra,'refcat_dec':matched_refcat_dec,
                                      'target_dist':target_dist,
                                      'gp1':{'mag':gp1_mag,'err':gp1_err},
                                      'rp1':{'mag':rp1_mag,'err':rp1_err},
                                      'ip1':{'mag':ip1_mag,'err':ip1_err},
                                      'zp1':{'mag':zp1_mag,'err':zp1_err},
                                      'g_sdss':{'mag':magnitudes['g_sdss']['mag'],'err':magnitudes['g_sdss']['err']},
                                      'r_sdss':{'mag':magnitudes['r_sdss']['mag'],'err':magnitudes['r_sdss']['err']},
                                      'i_sdss':{'mag':magnitudes['i_sdss']['mag'],'err':magnitudes['i_sdss']['err']},
                                      'z_sdss':{'mag':magnitudes['z_sdss']['mag'],'err':magnitudes['z_sdss']['err']},
                                      'B':{'mag':magnitudes['B']['mag'],'err':magnitudes['B']['err']},
                                      'V':{'mag':magnitudes['V']['mag'],'err':magnitudes['V']['err']},
                                      'R':{'mag':magnitudes['R']['mag'],'err':magnitudes['R']['err']},
                                      'I':{'mag':magnitudes['I']['mag'],'err':magnitudes['I']['err']}
                                     }
    os.remove(refcat_sources_filename_tmp)
    return refcat_photometry_dict

def generate_zeropoint_histogram(zeropoints,fits_filename,filename_suffix,param_dict):
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs_bgsub']
    mu,sigma = 0,0
    plot_filename = ext_filename_wcs_bgsub[:-5]+'.zeropoint_'+filename_suffix+'.pdf'
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
    plt.savefig(plot_filename,format='pdf',transparent=True,bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    output_log_entry(param_dict['log_filepath'],'Zero-point histogram generation and Gaussian fitting (mu={:.2f} mag; sigma={:.2f}) completed successfully.'.format(mu,sigma))
    return mu,sigma


def write_field_source_photometry_to_file(param_dict,fits_filename):
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs_bgsub']
    output_filename = ext_filename_wcs_bgsub[:-5] + '.field_source_photometry.txt'
    
    filter_name             = param_dict['fits_filenames'][fits_filename]['filter']
    field_sources           = param_dict['fits_filenames'][fits_filename]['field_sources']
    field_source_phot_table = param_dict['fits_filenames'][fits_filename]['field_source_phot_table']

    if len(field_sources) > 0 and len(field_source_phot_table) > 0:
        ra,dec                   = field_sources['ra'],field_sources['dec']
        xcoord,ycoord            = field_sources['xcentroid'],field_sources['ycentroid']
        flux,dflux               = field_source_phot_table['aperture_sum_bkgsub'],field_source_phot_table['aperture_sum_err']
        SNR                      = flux/dflux
        mag,magerr               = field_source_phot_table['mag'],field_source_phot_table['magerr']
                
        with open(output_filename,'w') as of:
            of.write('star_id    xcoord    ycoord             flux            dflux       SNR          RA         Dec     mag    err\n')
            num_ref_stars  = len(field_sources)
            for idx in range(0,num_ref_stars):
                if SNR[idx] > param_dict['daofind_sigma_threshold'] and mag[idx] != 99.999:
                    of.write('{:>7d}  {:8.3f}  {:8.3f}  {:15.3f}  {:15.3f}  {:8.1f}  {:10.6f}  {:10.6f}  {:6.3f}  {:5.3f}\n'.format(idx+1,
                               xcoord[idx],ycoord[idx],flux[idx],dflux[idx],SNR[idx],ra[idx],dec[idx],mag[idx],magerr[idx]))
                        
    return param_dict


def write_calibration_stars_to_file(param_dict,fits_filename):
    output_log_entry(param_dict['log_filepath'],'Writing calibration stars for {:s} to file...'.format(fits_filename))
    output_log_entry(param_dict['log_filepath'],'Writing calibration stars for absolute calibration of {:s} to file...'.format(fits_filename))
    ext_filename_wcs = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']
    
    if 'calib_star_table' in param_dict['fits_filenames'][fits_filename]:
        output_filename = ext_filename_wcs[:-5]+'.refcat_calib_stars.txt'
        calib_star_table = param_dict['fits_filenames'][fits_filename]['calib_star_table']
        #print(calib_star_table.columns)
    
        if calib_star_table != 0:
            with open(output_filename,'w') as of:
                of.write('star_id    xcoord    ycoord             flux            dflux       SNR          RA         Dec   refcat_RA  refcat_Dec  star_sep  gp1_mag  gp1_err  rp1_mag  rp1_err  ip1_mag  ip1_err  zp1_mag  zp1_err  g_sdss_mag  g_sdss_err  r_sdss_mag  r_sdss_err  i_sdss_mag  i_sdss_err  z_sdss_mag  z_sdss_err    B_mag    B_err    V_mag    V_err    R_mag    R_err    I_mag    I_err\n')
                num_calib_stars  = len(calib_star_table)
                for idx in range(0,num_calib_stars):
                    of.write('{:>7d}  {:8.3f}  {:8.3f}  {:15.3f}  {:15.3f}  {:8.1f}  {:10.6f}  {:10.6f}  {:10.6f}  {:10.6f}      {:4.1f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}      {:6.3f}       {:5.3f}      {:6.3f}       {:5.3f}      {:6.3f}       {:5.3f}      {:6.3f}       {:5.3f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}   {:6.3f}    {:5.3f}\n'.format(idx+1, \
                        calib_star_table[idx]['x'],calib_star_table[idx]['y'], \
                        calib_star_table[idx]['flux'],calib_star_table[idx]['dflux'],calib_star_table[idx]['SNR'], \
                        calib_star_table[idx]['ra'],calib_star_table[idx]['dec'], \
                        calib_star_table[idx]['refcat_ra'],calib_star_table[idx]['refcat_dec'],calib_star_table[idx]['target_dist'], \
                        calib_star_table[idx]['gp1_mag'],calib_star_table[idx]['gp1_err'], \
                        calib_star_table[idx]['rp1_mag'],calib_star_table[idx]['rp1_err'], \
                        calib_star_table[idx]['ip1_mag'],calib_star_table[idx]['ip1_err'], \
                        calib_star_table[idx]['zp1_mag'],calib_star_table[idx]['zp1_err'], \
                        calib_star_table[idx]['g_sdss_mag'],calib_star_table[idx]['g_sdss_err'], \
                        calib_star_table[idx]['r_sdss_mag'],calib_star_table[idx]['r_sdss_err'], \
                        calib_star_table[idx]['i_sdss_mag'],calib_star_table[idx]['i_sdss_err'], \
                        calib_star_table[idx]['z_sdss_mag'],calib_star_table[idx]['z_sdss_err'], \
                        calib_star_table[idx]['B_mag'],calib_star_table[idx]['B_err'], \
                        calib_star_table[idx]['V_mag'],calib_star_table[idx]['V_err'], \
                        calib_star_table[idx]['R_mag'],calib_star_table[idx]['R_err'], \
                        calib_star_table[idx]['I_mag'],calib_star_table[idx]['I_err']))

    return param_dict



def find_target_candidate_sources(param_dict,fits_filename):
    # Extracts field sources from image using DAOStarFinder
    # Based on code in https://photutils.readthedocs.io/en/stable/detection.html
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs'][:-5] + '_bgsub.fits'

    with fits.open(ext_filename_wcs_bgsub) as hdulist:
        data_bgsub      = hdulist[0].data
        ny,nx           = data_bgsub.shape
        mean,median,std = sigma_clipped_stats(data_bgsub,sigma=3.0)
        
    param_dict['fits_filenames'][fits_filename]['nx']     = nx
    param_dict['fits_filenames'][fits_filename]['ny']     = ny
    param_dict['fits_filenames'][fits_filename]['mean']   = mean
    param_dict['fits_filenames'][fits_filename]['median'] = median
    param_dict['fits_filenames'][fits_filename]['std']    = std
    
    daofind                  = DAOStarFinder(fwhm=param_dict['daofind_fwhm_target'],threshold=param_dict['daofind_sigma_threshold_target']*std)
    target_candidate_sources = daofind(data_bgsub)
    
    # if target source search fails, exit function
    if target_candidate_sources is None:
        num_target_candidate_sources = 0
        output_log_entry(param_dict['log_filepath'],'No target_candidate sources found')
        param_dict['fits_filenames'][fits_filename]['target_candidate_sources'] = target_candidate_sources
        return param_dict
    
    num_target_candidate_sources = len(target_candidate_sources['xcentroid'])
    header = fits.getheader(ext_filename_wcs_bgsub)
    w = WCS(header)
    target_candidate_sources['ra'],target_candidate_sources['dec'] = w.wcs_pix2world(target_candidate_sources['xcentroid'],target_candidate_sources['ycentroid'],1)
    target_candidate_sources['radec']                              = SkyCoord(target_candidate_sources['ra'],target_candidate_sources['dec'],unit='deg')
    target_candidate_sources['positions'] = np.transpose((target_candidate_sources['xcentroid'],target_candidate_sources['ycentroid']))
        
    # Do quick check of photometry
    field_source_phot_radius_pix = param_dict['field_source_phot_radius_pix']        
    positions              = target_candidate_sources['positions']
    aperture               = CircularAperture(positions,r=field_source_phot_radius_pix)
    data_bgsub,error,imwcs = extract_image_file_data(ext_filename_wcs_bgsub)
    sky_inner_r_pix        = param_dict['sky_inner_r_pix']
    sky_outer_r_pix        = param_dict['sky_outer_r_pix']
    annulus_aperture       = CircularAnnulus(positions,r_in=sky_inner_r_pix,r_out=sky_outer_r_pix)
    annulus_mask           = annulus_aperture.to_mask(method='center')
    bkg_median,bkg_stdev   = [],[]
    for mask in annulus_mask:
        annulus_data    = mask.multiply(data_bgsub)
        annulus_data_1d = annulus_data[mask.data > 0]
        _,median_sigclip,stdev_sigclip = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
        bkg_stdev.append(stdev_sigclip)
    bkg_median = np.array(bkg_median)
    bkg_stdev  = np.array(bkg_stdev)
    phot = aperture_photometry(data_bgsub,aperture,method='exact',error=error)
    phot['aperture_bkg']        = bkg_median*aperture.area
    phot['aperture_sum_bkgsub'] = phot['aperture_sum']-phot['aperture_bkg']
        
    # output source positions where fluxes are positive and flux > fluxerr
    target_candidate_sources_filename = ext_filename_wcs_bgsub[:-5]+'.target_candidate_sources.txt'
    with open(target_candidate_sources_filename,'w') as of:
        of.write('star_id     RA (deg)   Dec (deg)     xcoord    ycoord\n')
        for idx in range(len(target_candidate_sources['ra'])):
            if phot['aperture_sum_bkgsub'][idx] > 0 and phot['aperture_sum_bkgsub'][idx] > phot['aperture_sum_err'][idx]:
                of.write('  {:03d}    {:11.7f} {:11.7f}  {:9.3f} {:9.3f}\n'.format(idx+1,\
                            target_candidate_sources['ra'][idx],target_candidate_sources['dec'][idx],\
                            target_candidate_sources['xcentroid'][idx],target_candidate_sources['ycentroid'][idx]))
                    
    # Read filtered source positions back in
    target_candidate_sources = {}
    target_candidate_sources['source_idx'] = np.genfromtxt(target_candidate_sources_filename,skip_header=1,usecols=(0))
    target_candidate_sources['ra']         = np.genfromtxt(target_candidate_sources_filename,skip_header=1,usecols=(1))
    target_candidate_sources['dec']        = np.genfromtxt(target_candidate_sources_filename,skip_header=1,usecols=(2))
    target_candidate_sources['radec']      = SkyCoord(target_candidate_sources['ra'],target_candidate_sources['dec'],unit='deg')
    target_candidate_sources['xcentroid']  = np.genfromtxt(target_candidate_sources_filename,skip_header=1,usecols=(3))
    target_candidate_sources['ycentroid']  = np.genfromtxt(target_candidate_sources_filename,skip_header=1,usecols=(4))
    target_candidate_sources['positions']  = np.transpose((target_candidate_sources['xcentroid'],target_candidate_sources['ycentroid']))
    
    param_dict['fits_filenames'][fits_filename]['target_candidate_sources'] = target_candidate_sources
        
    output_log_entry(param_dict['log_filepath'],'{:d} target_candidate sources found'.format(num_target_candidate_sources))
    
    return param_dict


def plot_target_apertures(param_dict,fits_filename):
    # Plots target position and all apertures
    output_log_entry(param_dict['log_filepath'],'Plotting target photometry apertures for {:s}...'.format(fits_filename))
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs_bgsub']

    target_idx = param_dict['fits_filenames'][fits_filename]['target_idx']
    
    data_bgsub,error,imwcs = extract_image_file_data(ext_filename_wcs_bgsub)
    ny,nx = data_bgsub.shape
    
    target_phot_radii_pix  = param_dict['target_phot_radii_pix']
    sky_inner_r_pix        = param_dict['target_sky_inner_r_pix']
    sky_outer_r_pix        = param_dict['target_sky_outer_r_pix']

    target_xcentroid = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['xcentroid'][target_idx]
    target_ycentroid = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ycentroid'][target_idx]
    target_position  = np.transpose(([target_xcentroid],[target_ycentroid]))
    
    for target_phot_radius_pix in target_phot_radii_pix:
        aperture               = CircularAperture(target_position,r=target_phot_radius_pix)
        annulus_aperture       = CircularAnnulus(target_position,r_in=sky_inner_r_pix,r_out=sky_outer_r_pix)
    
        norm = simple_norm(data_bgsub,'sqrt',percent=99)
        plt.imshow(data_bgsub,norm=norm,interpolation='nearest')
        plt.xlim(target_xcentroid-sky_outer_r_pix*5,target_xcentroid+sky_outer_r_pix*5)
        plt.ylim(target_ycentroid-sky_outer_r_pix*5,target_ycentroid+sky_outer_r_pix*5)

        ap_patches = aperture.plot(color='white',lw=1,label='Photometry aperture')
        ann_patches = annulus_aperture.plot(color='red',lw=1,label='Background annulus')
        handles = (ap_patches[0],ann_patches[0])
        plt.legend(loc=(0.17,0.05),facecolor='#458989',labelcolor='white',handles=handles,prop={'weight':'bold','size':11})

        # Save to a File
        target_phot_radius_arcsec = target_phot_radius_pix * param_dict['instrument_params']['pixscale']
        plot_filepath = ext_filename_wcs_bgsub[:-5]+'.target_ap{:.1f}arcsec.pdf'.format(target_phot_radius_arcsec)
        plt.savefig(plot_filepath,format = 'pdf', transparent=True)
        plt.clf()
        plt.cla()
        plt.close()

        output_log_entry(param_dict['log_filepath'],'Target photometry apertures plotted ({:s}).'.format(plot_filepath))
    
    return param_dict


def find_target_source(param_dict,fits_filename):
    #param_dict = get_obj_ra_dec_horizons(param_dict,fits_filename)
    #param_dict = read_target_position_list(param_dict,fits_filename)
    param_dict = find_target_candidate_sources(param_dict,fits_filename)

    # Search for object using expected ephemeris position
    target_ra_horizons  = param_dict['fits_filenames'][fits_filename]['ephem']['ra_deg']
    target_dec_horizons = param_dict['fits_filenames'][fits_filename]['ephem']['dec_deg']
    target_x_horizons   = param_dict['fits_filenames'][fits_filename]['obj_x']
    target_y_horizons   = param_dict['fits_filenames'][fits_filename]['obj_y']
    output_log_entry(param_dict['log_filepath'],'Target expected at RA={:.5f},Dec={:.5f} (x={:.1f},y={:.1f})'.format(target_ra_horizons,target_dec_horizons,target_x_horizons,target_y_horizons))

    target_idx = -1
    for idx in range(len(param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ra'])):
        target_candidate_source_ra     = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ra'][idx]
        target_candidate_source_dec    = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['dec'][idx]
        x    = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['xcentroid'][idx]
        y    = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ycentroid'][idx]

        if abs(target_x_horizons-x) < 50 and abs(target_y_horizons-y) < 50:
            c1 = SkyCoord(target_candidate_source_ra,target_candidate_source_dec,unit='deg',frame='icrs')
            c2 = SkyCoord(target_ra_horizons,target_dec_horizons,unit='deg',frame='icrs')
            sep = c1.separation(c2)
            if sep.arcsecond < 5:
                target_idx = idx
                output_log_entry(param_dict['log_filepath'],'Target (RA={:.5f},Dec={:.5f}) found in target_candidate source list (idx={:d},RA={:.5f},Dec={:.5f})'.format(float(str(target_ra_horizons)),
                                        float(str(target_dec_horizons)),target_idx,float(str(target_candidate_source_ra)),float(str(target_candidate_source_dec))))
                output_log_entry(param_dict['log_filepath'],'({:.3f} arcsec from expected position)'.format(sep.arcsecond))
                
    if target_idx == -1:
        output_log_entry(param_dict['log_filepath'],'Target (RA={:.5f},Dec={:.5f}) not found in target_candidate source list'.format(float(str(target_ra_horizons)),float(str(target_dec_horizons))))
        
    param_dict['fits_filenames'][fits_filename]['target_idx'] = target_idx
    
    return param_dict


def measure_target_photometry(param_dict,fits_filename):
    target_idx = param_dict['fits_filenames'][fits_filename]['target_idx']
    #ext_filename_wcs = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']
    ext_filename_wcs_bgsub = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs'][:-5] + '_bgsub.fits'
    
    data_bgsub,error,imwcs = extract_image_file_data(ext_filename_wcs_bgsub)
    ny,nx = data_bgsub.shape
    
    pixscale                 = param_dict['instrument_params']['pixscale']
    #print(param_dict['target_phot_radii_arcsec'])
    target_phot_radii_arcsec = param_dict['target_phot_radii_arcsec']
    sky_inner_r_pix          = param_dict['target_sky_inner_r_arcsec'] / pixscale
    sky_outer_r_pix          = param_dict['target_sky_outer_r_arcsec'] / pixscale

    target_xcentroid = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['xcentroid'][target_idx]
    target_ycentroid = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ycentroid'][target_idx]
    target_position  = np.transpose(([target_xcentroid],[target_ycentroid]))
    
    param_dict['fits_filenames'][fits_filename]['target_phot_table'] = {}

    for target_phot_radius_arcsec in target_phot_radii_arcsec:
        target_phot_table = Table()
        target_phot_radius_pix = target_phot_radius_arcsec / pixscale
        #print('Aperture radius: {:.3f} --> {:.3f}'.format(target_phot_radius_arcsec,target_phot_radius_pix))
        aperture = CircularAperture(target_position,r=target_phot_radius_pix)
        
        # Compute image background properties from sky annulus
        annulus_aperture = CircularAnnulus(target_position,r_in=sky_inner_r_pix,r_out=sky_outer_r_pix)
        annulus_mask     = annulus_aperture.to_mask(method='center')
        bkg_median,bkg_stdev = [],[]
        for mask in annulus_mask:
            annulus_data    = mask.multiply(data_bgsub)
            annulus_data_1d = annulus_data[mask.data > 0]
            _,median_sigclip,stdev_sigclip = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            bkg_stdev.append(stdev_sigclip)
        bkg_median = np.array(bkg_median)
        bkg_stdev  = np.array(bkg_stdev)

        #print('bkg_stdev[0]:      {:.3f}'.format(bkg_stdev[0]))
        #print('median(bkg_stdev): {:.3f}'.format(statistics.median(bkg_stdev)))
        
        limit_mag_ps,limit_mag_sb = compute_detection_limits(param_dict,fits_filename,statistics.median(bkg_stdev),target_phot_radius_pix)
        
        # Run aperture photometry function
        phot = aperture_photometry(data_bgsub,aperture,method='exact',error=error)
        phot['aperture_bkg']        = bkg_median*aperture.area
        phot['aperture_sum_bkgsub'] = phot['aperture_sum']-phot['aperture_bkg']
    
        output_log_entry(param_dict['log_filepath'],'Using computed zero point from field stars: {:.3f}+/-{:.3f}'.format(param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag'],param_dict['fits_filenames'][fits_filename]['median_zeropoint_err']))
        zeropoint = ufloat(param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag'],param_dict['fits_filenames'][fits_filename]['median_zeropoint_err'])

        #print(zeropoint)
        
        exptime = param_dict['fits_filenames'][fits_filename]['exposure_time']
        phot['mag'],phot['magerr'] = fluxes2mags(phot['aperture_sum_bkgsub'],phot['aperture_sum_err'],zeropoint,exptime)
        
        #print(phot.colnames)
        #print(target_phot_table)
        
        # Transfer data from aperture photometry output table to phot_table and convert to mag
        apidx = 'ap{:.1f}arcsec'.format(target_phot_radius_arcsec)
        #target_phot_table.add_column(phot['aperture_sum'],name='aperture_sum')
        #target_phot_table.add_column(phot['aperture_sum_err'],name='aperture_sum_err')
        #target_phot_table.add_column(aperture.area,name='aperture_area')
        #target_phot_table.add_column(bkg_median,name='annulus_median')
        #target_phot_table.add_column(phot['aperture_bkg'],name='aperture_bkg')
        #target_phot_table.add_column((phot['aperture_sum_bkgsub']),name='aperture_sum_bkgsub')
        #target_phot_table.add_column((phot['mag']),name='mag')
        #target_phot_table.add_column((phot['magerr']),name='magerr')
        ##print(target_phot_table.colnames)
        #param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx] = target_phot_table
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx] = {}
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_sum']        = phot['aperture_sum'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_sum_err']    = phot['aperture_sum_err'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_area']       = aperture[0].area
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['annulus_median']      = bkg_median[0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['annulus_stdev']       = bkg_stdev[0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_bkg']        = phot['aperture_bkg'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_sum_bkgsub'] = phot['aperture_sum_bkgsub'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['mag']                 = phot['mag'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['magerr']              = phot['magerr'][0]
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['mag_corr']            = 0.0
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['magerr_corr']         = 0.0
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['limit_mag_ps']        = limit_mag_ps
        param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['limit_mag_sb']        = limit_mag_sb
        #print(target_phot_table.colnames)
        #param_dict['fits_filenames'][fits_filename]['target_phot_table'][apidx]['aperture_sum'] = phot['aperture_sum'][0]
        
        output_log_entry(param_dict['log_filepath'],'Target mag ({:.1f} arcsec/{:.2f} pix) = {:.3f}+/-{:.3f} mag'.format(target_phot_radius_arcsec,target_phot_radius_pix,phot['mag'][0],phot['magerr'][0]))
        output_log_entry(param_dict['log_filepath'],'Point source detection limit          = {:.3f} mag'.format(limit_mag_ps))
        output_log_entry(param_dict['log_filepath'],'Surface brightness detection limit    = {:.3f} mag/arcsec^2'.format(limit_mag_sb))
        output_log_entry(param_dict['log_filepath'],'Sky Standard Dev                      = {:.3f}'.format(statistics.median(bkg_stdev)))
        
    return param_dict



def compute_detection_limits(param_dict,fits_filename,sky_stddev,phot_radius_pix):
    limit_mag_ps,limit_mag_sb = 0,0
    exptime         = param_dict['fits_filenames'][fits_filename]['exposure_time']
    pixscale        = param_dict['instrument_params']['pixscale']
    zeropoint_mag   = param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag']
    ps_apert_noise = (math.pi*(phot_radius_pix)**2)**0.5 * sky_stddev  # using circular aperture with r=phot_radius_pix
    sb_apert_noise = ((1/pixscale)**2)**0.5 * sky_stddev                # using 1 arcsec x 1 arcsec aperture
    limit_flux_ps  = ps_apert_noise * 3
    limit_flux_sb  = sb_apert_noise * 3
    limit_mag_ps   = zeropoint_mag - 2.5*np.log10(limit_flux_ps/exptime)
    limit_mag_sb   = zeropoint_mag - 2.5*np.log10(limit_flux_sb/exptime)
    output_log_entry(param_dict['log_filepath'],'Point source detection limit:       {:.3f} mag'.format(limit_mag_ps))
    output_log_entry(param_dict['log_filepath'],'Surface brightness detection limit: {:.3f} mag/arcsec^2'.format(limit_mag_sb))
    return limit_mag_ps,limit_mag_sb


def clean_final_photometry_files(param_dict):
    os.chdir(param_dict['query_results_dirpath'])
    target_photometry_dirpath = param_dict['query_results_dirpath'] + 'target_photometry/'
    for filename in sorted(glob.glob('target_phot_results*')):
        os.rename(filename,target_photometry_dirpath+filename)
    return param_dict



def clean_output_files(param_dict,fits_filename):
    os.chdir(param_dict['query_results_dirpath'])

    orig_data_dirpath = param_dict['query_results_dirpath'] + 'data_files_orig/'
    if 'orig_data_dirpath' not in param_dict:
        param_dict['orig_data_dirpath'] = orig_data_dirpath

    ext_wcs_data_dirpath = param_dict['query_results_dirpath'] + 'data_files_ext_wcs/'
    if 'ext_wcs_data_dirpath' not in param_dict:
        param_dict['ext_wcs_data_dirpath'] = ext_wcs_data_dirpath

    ext_wcs_bgsub_data_dirpath = param_dict['query_results_dirpath'] + 'data_files_ext_wcs_bgsub/'
    if 'ext_wcs_bgsub_data_dirpath' not in param_dict:
        param_dict['ext_wcs_bgsub_data_dirpath'] = ext_wcs_bgsub_data_dirpath

    field_source_photometry_dirpath = param_dict['query_results_dirpath'] + 'field_source_photometry/'
    if 'field_source_photometry_dirpath' not in param_dict:
        param_dict['field_source_photometry_dirpath'] = field_source_photometry_dirpath

    target_photometry_dirpath = param_dict['query_results_dirpath'] + 'target_photometry/'
    if 'target_photometry_dirpath' not in param_dict:
        param_dict['target_photometry_dirpath'] = target_photometry_dirpath

    if not os.path.exists(orig_data_dirpath):               os.mkdir(orig_data_dirpath)
    if not os.path.exists(ext_wcs_data_dirpath):            os.mkdir(ext_wcs_data_dirpath)
    if not os.path.exists(ext_wcs_bgsub_data_dirpath):      os.mkdir(ext_wcs_bgsub_data_dirpath)
    if not os.path.exists(field_source_photometry_dirpath): os.mkdir(field_source_photometry_dirpath)
    if not os.path.exists(target_photometry_dirpath):       os.mkdir(target_photometry_dirpath)

    for filename in sorted(glob.glob(fits_filename[:-5]+'_??_wcs_bgsub.fits')):
        os.rename(filename,ext_wcs_bgsub_data_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'_??_wcs.fits')):
        os.rename(filename,ext_wcs_data_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'_??.fits')):
        os.rename(filename,ext_data_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename)):
        os.rename(filename,orig_data_dirpath+filename)

    for filename in sorted(glob.glob(fits_filename[:-5]+'*.field_source_photometry.txt')):
        os.rename(filename,field_source_photometry_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'*.field_sources.txt')):
        os.rename(filename,field_source_photometry_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'*.refcat_calib_stars.txt')):
        os.rename(filename,field_source_photometry_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'*.field_source_apertures.pdf')):
        os.rename(filename,field_source_photometry_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'*.zeropoint_histogram?.pdf')):
        os.rename(filename,field_source_photometry_dirpath+filename)

    for filename in sorted(glob.glob(fits_filename[:-5]+'*.target_candidate_sources.txt')):
        os.rename(filename,target_photometry_dirpath+filename)
    for filename in sorted(glob.glob(fits_filename[:-5]+'*arcsec.pdf')):
        os.rename(filename,target_photometry_dirpath+filename)

    os.chdir(orig_data_dirpath)
    for filename in sorted(glob.glob('*.fits')):
        fpack(filename,param_dict['config'])
    os.chdir(ext_wcs_data_dirpath)
    for filename in sorted(glob.glob('*.fits')):
        fpack(filename,param_dict['config'])
    os.chdir(ext_wcs_bgsub_data_dirpath)
    for filename in sorted(glob.glob('*.fits')):
        fpack(filename,param_dict['config'])
    
    os.chdir(param_dict['query_results_dirpath'])
    
    return param_dict



def write_target_photometry_to_file(param_dict,fits_filename):
    target_phot_radii_arcsec = param_dict['target_phot_radii_arcsec']
    filters = ['g_sdss','r_sdss','i_sdss','z_sdss','B','V','R','I']
    for target_phot_radius_arcsec in target_phot_radii_arcsec:
        output_filename = 'target_phot_results.ap{:04.1f}arcsec.txt'.format(target_phot_radius_arcsec)
        target = param_dict['target_name']
        telescope = param_dict['telinst']

        with open(output_filename,'a') as of:
            filter_name = param_dict['fits_filenames'][fits_filename]['filter']
            target_ext_id = param_dict['fits_filenames'][fits_filename]['target_ext_id']
            if filter_name in filters and target_ext_id != -1:
                time_mid       = param_dict['fits_filenames'][fits_filename]['time_mid']
                jdmid          = time_mid.jd
                utdate         = time_mid.iso[:10]
                uttimemid      = time_mid.iso[11:21]
                uttimemidhr    = int(uttimemid[:2]) + int(uttimemid[3:5])/60 + float(uttimemid[6:])/3600
                exptime        = param_dict['fits_filenames'][fits_filename]['exposure_time']
                filter_name    = param_dict['fits_filenames'][fits_filename]['filter']
                airmass        = param_dict['fits_filenames'][fits_filename]['ephem']['airmass']
                ephemra        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_deg']
                ephemdec       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_deg']
                ra_rate        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_rate']
                dec_rate       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_rate']
                tot_rate       = param_dict['fits_filenames'][fits_filename]['ephem']['tot_rate']
                expt_trail_pa  = param_dict['fits_filenames'][fits_filename]['ephem']['trail_pa_expected']
                expt_trail_len = param_dict['fits_filenames'][fits_filename]['ephem']['trail_len_expected']
                rdist          = param_dict['fits_filenames'][fits_filename]['ephem']['rdist']
                ddist          = param_dict['fits_filenames'][fits_filename]['ephem']['ddist']
                phsang         = param_dict['fits_filenames'][fits_filename]['ephem']['phase_angle']
                pa_as          = param_dict['fits_filenames'][fits_filename]['ephem']['pa_as']
                pa_nhv         = param_dict['fits_filenames'][fits_filename]['ephem']['pa_nhv']
                orbplang       = param_dict['fits_filenames'][fits_filename]['ephem']['orb_pl_angle']
                trueanom       = param_dict['fits_filenames'][fits_filename]['ephem']['true_anomaly']
                v_mag_str      = param_dict['fits_filenames'][fits_filename]['ephem']['v_mag_str']
                n_mag_str      = param_dict['fits_filenames'][fits_filename]['ephem']['n_mag_str']
                t_mag_str      = param_dict['fits_filenames'][fits_filename]['ephem']['t_mag_str']
                ra_3sig        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_3sigma']
                dec_3sig       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_3sigma']
            
                if 'target_idx' in param_dict['fits_filenames'][fits_filename]:
                    target_idx  = param_dict['fits_filenames'][fits_filename]['target_idx']
                else:
                    target_idx = -1
                
                if target_idx == -1:
                    xcoord,ycoord,source_ra,source_dec = 0,0,0,0
                    flux,dflux,snr,mag,magerr = 0,0,0,0,0
                else:
                    xcoord      = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['xcentroid'][target_idx]
                    ycoord      = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ycentroid'][target_idx]
                    source_ra   = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ra'][target_idx]
                    source_dec  = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['dec'][target_idx]
                    flux        = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['aperture_sum_bkgsub']
                    dflux       = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['aperture_sum_err']
                    snr         = flux/dflux
                    mag         = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['mag']
                    magerr      = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['magerr']

                if 'target_phot_table' in param_dict['fits_filenames'][fits_filename]:
                    limit_mag_ps = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['limit_mag_ps']
                    limit_mag_sb = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['limit_mag_sb']
                else:
                    limit_mag_ps = 0.0
                    limit_mag_sb = 0.0

                if 'calib_star_table' in param_dict['fits_filenames'][fits_filename]:
                    nrefstars = len(param_dict['fits_filenames'][fits_filename]['calib_star_table'])
                else:
                    nrefstars = 0
                if 'median_zeropoint_mag' in param_dict['fits_filenames'][fits_filename]:
                    zeropt = param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag']
                    zpterr = param_dict['fits_filenames'][fits_filename]['median_zeropoint_err']
                else:
                    if 'zeropoint' in param_dict['fits_filenames'][fits_filename]:
                        zeropt = param_dict['fits_filenames'][fits_filename]['zeropoint']
                    else:
                        zeropt = 0.0
                    zpterr = 0.0

                of.write('{:<20s}  {:<20s}  {:<50s}  {:14.6f}  '.format(target,telescope,fits_filename,jdmid))
                of.write('{:10s}  {:12s}  {:11.6f}  {:7.1f}  '.format(utdate,uttimemid,uttimemidhr,exptime))
                of.write('{:<6s}  {:7.4f}  {:11.7f}  {:11.7f}  '.format(filter_name,airmass,ephemra,ephemdec))
                of.write('{:12.4f}  {:12.4f}  {:13.4f}  {:13.1f}  {:14.3f}  {:9.4f}  {:8.4f}  '.format(ra_rate,dec_rate,tot_rate,expt_trail_pa,expt_trail_len,rdist,ddist))
                of.write('{:6.2f}  {:7.3f}  {:7.3f}  {:8.3f}  '.format(phsang,pa_as,pa_nhv,orbplang))
                of.write('{:8.3f}  {:s}  {:s}  {:s}  {:9.4f}  {:10.4f}  '.format(trueanom,v_mag_str,n_mag_str,t_mag_str,ra_3sig,dec_3sig))
                of.write('{:9.3f}  {:9.3f}  {:15.3f}  {:15.3f}  '.format(xcoord,ycoord,flux,dflux))
                of.write('{:9.1f}  {:11.7f}  {:11.7f}  {:6.3f}  {:6.3f}  '.format(snr,source_ra,source_dec,zeropt,zpterr))
                of.write('{:>9d}  {:5.1f}  {:6.3f}  {:6.3f}  '.format(nrefstars,target_phot_radius_arcsec,mag,magerr))
                of.write('{:8.3f}  {:8.3f}  '.format(limit_mag_ps,limit_mag_sb))
                of.write('\n')
                    
    output_log_entry(param_dict['log_filepath'],'Photometry results written to {:s}.'.format(output_filename))
    return param_dict

def write_target_photometry_to_file_orig(param_dict):
    target_phot_radii_arcsec = param_dict['target_phot_radii_arcsec']
    filters = ['g_sdss','r_sdss','i_sdss','z_sdss','B','V','R','I']
    for target_phot_radius_arcsec in target_phot_radii_arcsec:
        output_filename = 'target_phot_results.ap{:04.1f}arcsec.txt'.format(target_phot_radius_arcsec)
        target = param_dict['target_name']
        telescope = param_dict['telinst']
        with open(output_filename,'w') as of:
            of.write('target                telescope             filename                                                     jdmid  ')
            of.write('UTdate      UTtimemid     UTtimemidHr  exptime  ')
            of.write('filter  airmass      EphemRA     EphemDec  ')
            of.write('     RA_rate      Dec_rate       tot_rate  expt_trail_pa  expt_trail_len  heliodist   geodist  ')
            of.write('phsang    PA_AS   PA_NHV  orbplang  ')
            of.write('trueanom  RA_3sigma  Dec_3sigma  ')
            of.write('   xcoord     ycoord             flux            dflux   ')
            of.write('     SNR           RA          Dec  zeropt  zpterr  ')
            of.write('nrefstars  apert     mag  magerr  limit_ps  limit_sb')
            of.write('\n')

            for fits_filename in param_dict['fits_filenames']:
                filter_name = param_dict['fits_filenames'][fits_filename]['filter']
                target_ext_id = param_dict['fits_filenames'][fits_filename]['target_ext_id']
                if filter_name in filters and target_ext_id != -1:
                    time_mid       = param_dict['fits_filenames'][fits_filename]['time_mid']
                    jdmid          = time_mid.jd
                    utdate         = time_mid.iso[:10]
                    uttimemid      = time_mid.iso[11:21]
                    uttimemidhr    = int(uttimemid[:2]) + int(uttimemid[3:5])/60 + float(uttimemid[6:])/3600
                    exptime        = param_dict['fits_filenames'][fits_filename]['exposure_time']
                    filter_name    = param_dict['fits_filenames'][fits_filename]['filter']
                    airmass        = param_dict['fits_filenames'][fits_filename]['ephem']['airmass']
                    ephemra        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_deg']
                    ephemdec       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_deg']
                    ra_rate        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_rate']
                    dec_rate       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_rate']
                    tot_rate       = param_dict['fits_filenames'][fits_filename]['ephem']['tot_rate']
                    expt_trail_pa  = param_dict['fits_filenames'][fits_filename]['ephem']['trail_pa_expected']
                    expt_trail_len = param_dict['fits_filenames'][fits_filename]['ephem']['trail_len_expected']
                    rdist          = param_dict['fits_filenames'][fits_filename]['ephem']['rdist']
                    ddist          = param_dict['fits_filenames'][fits_filename]['ephem']['ddist']
                    phsang         = param_dict['fits_filenames'][fits_filename]['ephem']['phase_angle']
                    pa_as          = param_dict['fits_filenames'][fits_filename]['ephem']['pa_as']
                    pa_nhv         = param_dict['fits_filenames'][fits_filename]['ephem']['pa_nhv']
                    orbplang       = param_dict['fits_filenames'][fits_filename]['ephem']['orb_pl_angle']
                    trueanom       = param_dict['fits_filenames'][fits_filename]['ephem']['true_anomaly']
                    ra_3sig        = param_dict['fits_filenames'][fits_filename]['ephem']['ra_3sigma']
                    dec_3sig       = param_dict['fits_filenames'][fits_filename]['ephem']['dec_3sigma']
                
                    if 'target_idx' in param_dict['fits_filenames'][fits_filename]:
                        target_idx  = param_dict['fits_filenames'][fits_filename]['target_idx']
                    else:
                        target_idx = -1
                    
                    if target_idx == -1:
                        xcoord,ycoord,source_ra,source_dec = 0,0,0,0
                        flux,dflux,snr,mag,magerr = 0,0,0,0,0
                    else:
                        xcoord      = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['xcentroid'][target_idx]
                        ycoord      = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ycentroid'][target_idx]
                        source_ra   = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['ra'][target_idx]
                        source_dec  = param_dict['fits_filenames'][fits_filename]['target_candidate_sources']['dec'][target_idx]
                        flux        = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['aperture_sum_bkgsub']
                        dflux       = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['aperture_sum_err']
                        snr         = flux/dflux
                        mag         = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['mag']
                        magerr      = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['magerr']

                    if 'target_phot_table' in param_dict['fits_filenames'][fits_filename]:
                        limit_mag_ps = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['limit_mag_ps']
                        limit_mag_sb = param_dict['fits_filenames'][fits_filename]['target_phot_table']['ap{:.1f}arcsec'.format(target_phot_radius_arcsec)]['limit_mag_sb']
                    else:
                        limit_mag_ps = 0.0
                        limit_mag_sb = 0.0

                    if 'calib_star_table' in param_dict['fits_filenames'][fits_filename]:
                        nrefstars = len(param_dict['fits_filenames'][fits_filename]['calib_star_table'])
                    else:
                        nrefstars = 0
                    if 'median_zeropoint_mag' in param_dict['fits_filenames'][fits_filename]:
                        zeropt = param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag']
                        zpterr = param_dict['fits_filenames'][fits_filename]['median_zeropoint_err']
                    else:
                        if 'zeropoint' in param_dict['fits_filenames'][fits_filename]:
                            zeropt = param_dict['fits_filenames'][fits_filename]['zeropoint']
                        else:
                            zeropt = 0.0
                        zpterr = 0.0

                    of.write('{:<20s}  {:<20s}  {:<50s}  {:14.6f}  '.format(target,telescope,fits_filename,jdmid))
                    of.write('{:10s}  {:12s}  {:11.6f}  {:7.1f}  '.format(utdate,uttimemid,uttimemidhr,exptime))
                    of.write('{:<6s}  {:7.4f}  {:11.7f}  {:11.7f}  '.format(filter_name,airmass,ephemra,ephemdec))
                    of.write('{:12.4f}  {:12.4f}  {:13.4f}  {:13.1f}  {:14.3f}  {:9.4f}  {:8.4f}  '.format(ra_rate,dec_rate,tot_rate,expt_trail_pa,expt_trail_len,rdist,ddist))
                    of.write('{:6.2f}  {:7.3f}  {:7.3f}  {:8.3f}  '.format(phsang,pa_as,pa_nhv,orbplang))
                    of.write('{:8.3f}  {:9.4f}  {:10.4f}  '.format(trueanom,ra_3sig,dec_3sig))
                    of.write('{:9.3f}  {:9.3f}  {:15.3f}  {:15.3f}  '.format(xcoord,ycoord,flux,dflux))
                    of.write('{:9.1f}  {:11.7f}  {:11.7f}  {:6.3f}  {:6.3f}  '.format(snr,source_ra,source_dec,zeropt,zpterr))
                    of.write('{:>9d}  {:5.1f}  {:6.3f}  {:6.3f}  '.format(nrefstars,target_phot_radius_arcsec,mag,magerr))
                    of.write('{:8.3f}  {:8.3f}  '.format(limit_mag_ps,limit_mag_sb))
                    of.write('\n')
                    
    output_log_entry(param_dict['log_filepath'],'Photometry results written to {:s}.'.format(output_filename))
    return param_dict



###################### OLD FUNCTIONS ######################

#def perform_photometry(param_dict):
#    for fits_filename in param_dict['fits_filenames']:
#        if param_dict['fits_filenames'][fits_filename]['target_ext_id'] != -1:
#            ext_filename_wcs = fits_filename[:-5] + '_{:02d}_wcs.fits'.format(ext_id)
#            ext_filename_wcs_bgsub = ext_filename_wcs[:-5] + '_bgsub.fits'
#            with fits.open(ext_filename_wcs) as hdu_data:
#                data          = hdu_data[0].data
#                hdr           = hdu_data[0].header
#                sigma_clip    = SigmaClip(sigma=3.0)
#                bkg_estimator = MedianBackground()
#                bkg = Background2D(data,(50,50),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
#                data_bgsub    = data - bkg.background
#                fits.writeto(ext_filename_wcs_bgsub,data_bgsub,hdr,overwrite=True,checksum=True)
#            psf_model = IntegratedGaussianPRF(flux=1,sigma=2.7/2.35)
#            fit_shape = (5,5)
#            finder = DAOStarFinder(6.0,2.0)
#            psfphot = PSFPhotometry(psf_model,fit_shape,finder=finder,aperture_radius=4)
#            phot = psfphot(data_bgsub,error=bkg.background_rms_median)
#
#    return None

# def compute_mags_from_zpoint(param_dict,fits_filename,target_flux,target_fluxerr):
#     # Compute magnitude of source from zero point
#     exptime               = param_dict['fits_filenames'][fits_filename]['exposure_time']
#     median_zeropoint_mag  = param_dict['fits_filenames'][fits_filename]['median_zeropoint_mag']
#     median_zeropoint_err  = param_dict['fits_filenames'][fits_filename]['median_zeropoint_err']
#     target_mag,target_err = 0,0
#     target_flux_ufloat    = ufloat(target_flux,target_fluxerr)
#     zpoint_mag_ufloat     = ufloat(zeropoint_mag,zeropoint_magerr)
#     target_mag_ufloat     = zpoint_mag_ufloat - 2.5 * unumpy.log10(target_flux_ufloat/exptime)
#     target_mag            = target_mag_ufloat.nominal_value
#     target_err            = target_mag_ufloat.std_dev
#     return target_mag,target_err

# def check_if_radec_match(ra1,dec1,ra2,dec2,threshold_arcsec):
#     radec_matched = False
#     test_threshold = threshold_arcsec/3600*3
#     if abs(dec1-dec2) < test_threshold and (abs(ra1-ra2) < test_threshold or abs(abs(ra1-ra2)-360) < test_threshold):
#         c1  = SkyCoord(ra1,dec1,unit='deg',frame='icrs')
#         c2  = SkyCoord(ra2,dec2,unit='deg',frame='icrs')
#         sep = c1.separation(c2)
#         if sep.arcsecond < threshold_arcsec:
#             radec_matched = True
#     return radec_matched

##### REFCAT DATA RETRIEVAL FUNCTION #####

# def retrieve_refcat_sources(param_dict,fits_filename):
#     # Retrieve refcat sources
#     ext_filename_wcs = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']
#     
#     output_log_entry(param_dict['log_filepath'],'Retrieving refcat sources for {:s}...'.format(ext_filename_wcs))
#     refcat_sources_filename_tmp = ext_filename_wcs[:-5]+'.refcat_sources_temp.txt'
#     refcat_sources_filename     = ext_filename_wcs[:-5]+'.refcat_sources.txt'
# 
#     # Determine search range
#     nx = param_dict['fits_filenames'][fits_filename]['npix_x']
#     ny = param_dict['fits_filenames'][fits_filename]['npix_y']
#     x_mid,y_mid = int(nx/2),int(ny/2)
#     header = fits.getheader(ext_filename_wcs)
#     w = WCS(header)
#     ra_mid,dec_mid = w.wcs_pix2world(x_mid,y_mid,1)
#     ra_0,dec_0     = w.wcs_pix2world(1,1,1)
#     c1 = SkyCoord(ra_mid,dec_mid,unit='deg',frame='icrs')
#     c2 = SkyCoord(ra_0,dec_0,unit='deg',frame='icrs')
#     sep = c1.separation(c2)
#     radius_deg = sep.degree
#     output_log_entry(param_dict['log_filepath'],'Search region: RA={:.5f} deg, Dec={:.5f} deg, radius={:.1f} arcmin'.format(ra_mid,dec_mid,radius_deg*60))
#     
#     with open(refcat_sources_filename_tmp,'w') as outf:
#         cmd_ra     = '{:f}'.format(ra_mid)
#         cmd_dec    = '{:f}'.format(dec_mid)
#         cmd_radius = '{:f}'.format(radius_deg)
#         cmd_refcat = param_dict['config']['cmd_refcat']
#         directories = '{:s},{:s},{:s},{:s},{:s}'.format(param_dict['config']['refcat_dir_00_m_16'],
#                         param_dict['config']['refcat_dir_16_m_17'],param_dict['config']['refcat_dir_17_m_18'],
#                         param_dict['config']['refcat_dir_18_m_19'],param_dict['config']['refcat_dir_19_m_20'])
#         #print(directories)
#         var_fields = 'RA,Dec,g,dg,r,dr,i,di,z,dz'
#         cmd_dir,cmd_rad,cmd_var,cmd_nohdr,cmd_bin = '-dir','-rad','-var','-nohdr','-bin'
#         cmd = [cmd_refcat,cmd_ra,cmd_dec,cmd_dir,directories,cmd_bin,cmd_rad,cmd_radius,cmd_var,var_fields,cmd_nohdr]
#         #print(cmd)
#         subprocess.call(cmd,stdout=outf)
#     with open(refcat_sources_filename,'w') as of, open(refcat_sources_filename_tmp,'r') as infile:
#         of.write(' objRA       objDec      gmag   gerr  rmag   rerr  imag   ierr  zmag   zerr\n')
#         for line in infile:
#             of.write(line)
#     os.remove(refcat_sources_filename_tmp)
#     output_log_entry(param_dict['log_filepath'],'refcat catalog sources for {:s} written to {:s}.'.format(ext_filename_wcs,refcat_sources_filename))
#         
#     # parse local file into astropy.table object
#     objRA  = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(0))
#     objDec = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(1))
#     gmag   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(2))
#     gerr   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(3))
#     rmag   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(4))
#     rerr   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(5))
#     imag   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(6))
#     ierr   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(7))
#     zmag   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(8))
#     zerr   = np.genfromtxt(refcat_sources_filename,skip_header=1,usecols=(9))
#     refcat_table = Table([objRA,objDec,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr], \
#         names=('right_ascension','declination','gp1_mag','gp1_err','rp1_mag','rp1_err','ip1_mag','ip1_err','zp1_mag','zp1_err'), \
#         dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
#     num_refcat_sources = len(refcat_table)
#     if num_refcat_sources > 0:
#         idx = 0
#         # Remove sources that are likely to be saturated or have dummy values (i.e., -999)
#         while idx < len(refcat_table):
#             if refcat_table[idx]['gp1_mag'] < 10 or refcat_table[idx]['rp1_mag'] < 10 or refcat_table[idx]['ip1_mag'] < 10 or refcat_table[idx]['zp1_mag'] < 10:
#                 refcat_table.remove_row(idx)
#             else:
#                 idx += 1
#         num_refcat_sources = len(refcat_table)
#     if num_refcat_sources > 0:
#         output_log_entry(param_dict['log_filepath'],'{:d} refcat catalog sources retrieved for {:s}.'.format(num_refcat_sources,ext_filename_wcs))
#     else:
#         output_log_entry(param_dict['log_filepath'],'No refcat calibration sources found for {:s}'.format(ext_filename_wcs))
#     
#     param_dict['fits_filenames'][fits_filename]['refcat_table'] = refcat_table
#         
#     return param_dict
#
#
# def match_refcat_sources(param_dict,fits_filename):
#     # Find matches for sources extracted by DAOStarFinder in refcat source catalog
#     output_log_entry(param_dict['log_filepath'],'Finding matches for sources extracted by DAOStarFinder in refcat catalog...')
#     ext_filename_wcs   = param_dict['fits_filenames'][fits_filename]['ext_filename_wcs']
#     
#     refcat_table       = param_dict['fits_filenames'][fits_filename]['refcat_table']
#     field_source_table = param_dict['fits_filenames'][fits_filename]['field_source_table']
#     calib_star_table   = Table(names=('x','y','ra','dec','flux','dflux','SNR','refcat_ra','refcat_dec',\
#                                       'gp1_mag','gp1_err','rp1_mag','rp1_err',\
#                                       'ip1_mag','ip1_err','zp1_mag','zp1_err',\
#                                       'g_sdss_mag','g_sdss_err','r_sdss_mag','r_sdss_err',\
#                                       'i_sdss_mag','i_sdss_err','z_sdss_mag','z_sdss_err',\
#                                       'B_mag','B_err','V_mag','V_err','R_mag','R_err','I_mag','I_err',\
#                                       'target_dist','zero_point','zpoint_err'))
#     num_field_sources = len(field_source_table)
#     num_refcat_stars  = len(refcat_table)
#     
#     output_log_entry(param_dict['log_filepath'],'Matching field sources to refcat catalog sources...')
#     stars_matched = 0
#     #print(field_source_table.columns)
#     #print(calib_star_table.columns)
#     for idx1 in range(0,num_field_sources):
#         idx2,star_matched = 0,0
#         while idx2 < num_refcat_stars and star_matched == 0:
#             ra1_deg  = param_dict['fits_filenames'][fits_filename]['field_sources']['ra'][idx1]
#             dec1_deg = param_dict['fits_filenames'][fits_filename]['field_sources']['dec'][idx1]
#             ra2_deg  = refcat_table[idx2]['right_ascension']
#             dec2_deg = refcat_table[idx2]['declination']
#             
#             if abs(dec1_deg-dec2_deg) < 0.005 and (abs(ra1_deg-ra2_deg)<0.005 or abs(abs(ra1_deg-ra2_deg)-360)<0.005):
#                 c1 = param_dict['fits_filenames'][fits_filename]['field_sources']['radec'][idx1]
#                 c2 = SkyCoord(refcat_table[idx2]['right_ascension'],refcat_table[idx2]['declination'],unit='deg',frame='icrs')
#                 sep = c1.separation(c2)
#                 
#                 if sep.arcsecond < 2.6:  # If target matches within 2.5 arcsec
#                     star_matched = 1
#                     copy_fieldsource_row_to_calibstar_table(field_source_table,calib_star_table,idx1)
#                     calib_star_table[stars_matched]['refcat_ra']   = float(refcat_table[idx2]['right_ascension'])
#                     calib_star_table[stars_matched]['refcat_dec']  = float(refcat_table[idx2]['declination'])
#                     calib_star_table[stars_matched]['target_dist'] = sep.arcsecond
#                     calib_star_table[stars_matched]['gp1_mag']     = refcat_table[idx2]['gp1_mag']
#                     calib_star_table[stars_matched]['gp1_err']     = refcat_table[idx2]['gp1_err']
#                     calib_star_table[stars_matched]['rp1_mag']     = refcat_table[idx2]['rp1_mag']
#                     calib_star_table[stars_matched]['rp1_err']     = refcat_table[idx2]['rp1_err']
#                     calib_star_table[stars_matched]['ip1_mag']     = refcat_table[idx2]['ip1_mag']
#                     calib_star_table[stars_matched]['ip1_err']     = refcat_table[idx2]['ip1_err']
#                     calib_star_table[stars_matched]['zp1_mag']     = refcat_table[idx2]['zp1_mag']
#                     calib_star_table[stars_matched]['zp1_err']     = refcat_table[idx2]['zp1_err']
#                     stars_matched += 1
#                     output_log_entry(param_dict['log_filepath'],'Field source at {:.6f},{:.6f} --> refcat source at {:.6f},{:.6f}'.format(ra1_deg,dec1_deg,ra2_deg,dec2_deg))
#                 else:
#                     idx2 += 1
#             else:
#                 idx2 += 1
#     output_log_entry(param_dict['log_filepath'],'{:d} field sources matched to refcat catalog sources for {:s}.'.format(stars_matched,ext_filename_wcs))
#     if stars_matched == 0:
#         output_log_entry(param_dict['log_filepath'],'No refcat sources matched to field sources for {:s}.'.format(ext_filename_wcs))
#         calib_star_table = 0
#     elif stars_matched < 3:
#         output_log_entry(param_dict['log_filepath'],'Insufficient refcat sources (n<3) matched to field sources for {:s}.'.format(ext_filename_wcs))
#         calib_star_table = 0
#         
#     param_dict['fits_filenames'][fits_filename]['calib_star_table'] = calib_star_table
#     
#     return param_dict
# 
# def copy_fieldsource_row_to_calibstar_table(field_source_table,calib_star_table,target_idx):
#     # Copy source from refcat source list to calibration star list
#     calib_star_table.add_row((field_source_table[target_idx]['x'],field_source_table[target_idx]['y'], \
#         field_source_table[target_idx]['ra'],field_source_table[target_idx]['dec'], \
#         field_source_table[target_idx]['flux'],field_source_table[target_idx]['dflux'], \
#         field_source_table[target_idx]['SNR'], \
#         field_source_table[target_idx]['refcat_ra'],field_source_table[target_idx]['refcat_dec'], \
#         field_source_table[target_idx]['gp1_mag'],field_source_table[target_idx]['gp1_err'], \
#         field_source_table[target_idx]['rp1_mag'],field_source_table[target_idx]['rp1_err'], \
#         field_source_table[target_idx]['ip1_mag'],field_source_table[target_idx]['ip1_err'], \
#         field_source_table[target_idx]['zp1_mag'],field_source_table[target_idx]['zp1_err'], \
#         field_source_table[target_idx]['g_sdss_mag'],field_source_table[target_idx]['g_sdss_err'], \
#         field_source_table[target_idx]['r_sdss_mag'],field_source_table[target_idx]['r_sdss_err'], \
#         field_source_table[target_idx]['i_sdss_mag'],field_source_table[target_idx]['i_sdss_err'], \
#         field_source_table[target_idx]['z_sdss_mag'],field_source_table[target_idx]['z_sdss_err'], \
#         field_source_table[target_idx]['B_mag'],field_source_table[target_idx]['B_err'], \
#         field_source_table[target_idx]['V_mag'],field_source_table[target_idx]['V_err'], \
#         field_source_table[target_idx]['R_mag'],field_source_table[target_idx]['R_err'], \
#         field_source_table[target_idx]['I_mag'],field_source_table[target_idx]['I_err'], \
#         field_source_table[target_idx]['target_dist'], \
#         field_source_table[target_idx]['zero_point'],field_source_table[target_idx]['zpoint_err']))
#     return calib_star_table
#
# ##### FUNCTION DEFINITIONS -- BASIC FILE FUNCTIONS #####
# 
# def create_directory(path,path_logfile,path_errorfile):
#     if not os.path.exists(path):
#         try:
#             os.mkdir(path)
#         except OSError:
#             output_error_log_entry(path_logfile,path_errorfile,'Creation of directory {:s} failed.'.format(path))
#         else:
#             output_log_entry(path_logfile,'Directory {:s} successfully created.'.format(path))
#     return None
# 
# def rm_file(path_file):
#     try:
#         cmd_rm = '/bin/rm'
#         cmd = [cmd_rm,path_file]
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: rm_file()')
#         print(e)
#     return None
# 
# def rename_file(old_filename,new_filename):
#     try:
#         cmd_mv = '/bin/mv'
#         cmd = [cmd_mv,old_filename,new_filename]
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: rename_file')
#         print(e)
#     return None
# 
# def move_file(old_filename,new_filename):
#     try:
#         cmd_mv = '/bin/mv'
#         cmd = [cmd_mv,old_filename,new_filename]
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: move_file')
#         print(e)
#     return None
# 
# def copy_file(old_filename,new_filepath):
#     try:
#         cmd_cp = '/bin/cp'
#         cmd = [cmd_cp,old_filename,new_filepath]
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: copy_file')
#         print(e)
#     return None
# 
# def decompress_file_bzip2(path_file):
#     try:
#         bzip1,bzip2 = 'bzip2','-d'
#         cmd = [bzip1,bzip2,path_file]
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: decompress_file_bzip2')
#         print(e)
#     return None
# 
# def compress_file_gzip(path_file):
#     try:
#         cmd_gzip = 'gzip'
#         cmd = [cmd_gzip,path_file]
#         if os.path.exists(path_file+'.gz'):
#             os.rename(path_file+'.gz',path_file+'.backup.gz')
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: compress_file_gzip')
#         print(e)
#     return None
# 
# def decompress_file_gzip(path_file):
#     try:
#         cmd_gunzip = 'gunzip'
#         cmd = [cmd_gunzip,path_file]
#         if os.path.exists(path_file[:-3]) and path_file[-3:] == '.gz':
#             os.rename(path_file[:-3],path_file[:-3]+'.backup')
#         process = subprocess.call(cmd)
#     except Exception as e:
#         print('Function failed: decompress_file_gzip')
#         print(e)
#     return None
# 
# def sky2pix_fzfile(fits_file,ra,dec):
#     xpix,ypix = -1,-1
#     try:
#         cmd_sky2pix = ''
#         if os.path.isfile('/usr/local/bin/sky2pix'):
#             cmd_sky2pix = '/usr/local/bin/sky2pix'
#         elif os.path.isfile('/users/hhsieh/dropbox/hhsieh_db/astro_db/tools/wcseval/sky2pix'):
#             cmd_sky2pix = '/users/hhsieh/dropbox/hhsieh_db/astro_db/tools/wcseval/sky2pix'
#         elif os.path.isfile('/atlas/bin/sky2pix'):
#             cmd_sky2pix = '/atlas/bin/sky2pix'
#         if cmd_sky2pix != '':
#             ra_str  = '{:.10f}'.format(ra)
#             dec_str = '{:.10f}'.format(dec)
#             cmd = [cmd_sky2pix,fits_file,ra_str,dec_str]
#             process = subprocess.Popen(cmd,stdout=subprocess.PIPE)
#             out,err = process.communicate()
#             output = out.split()
#             xpix = float(output[0])
#             ypix = float(output[1])
#         else:
#             print('Function failed: sky2pix_fzfile (command not found)')
#     except Exception as e:
#         print('Function failed: sky2pix_fzfile')
#         print(e)
#     return xpix,ypix
# 
# 
# def compute_angular_dist_radec(ra1,dec1,ra2,dec2):
#     # Compute angular distance in degrees from two sets of RA,Dec coordinates, also in degrees
#     source_distance_deg = -1
#     try:
#         c1 = SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
#         c2 = SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
#         source_distance_deg = c1.separation(c2).degree
#     except Exception as e:
#         print('Function failed: compute_angular_dist()')
#         print(e)
#     return source_distance_deg
# 
# 
# def get_pixel_coords_from_radec(fits_file,ra,dec):
#     px,py = [-1,-1],[-1,-1]
#     try:
#         header = fits.getheader(fits_file)
#         w = WCS(header)
#         px,py = w.wcs_world2pix(ra,dec,1)
#     except Exception as e:
#         print('Function failed: get_pixel_coords_from_radec()')
#         print(e)
#     return px,py
# 
# 
# def get_radec_from_pixel_coords(fits_file,xcoord,ycoord):
#     ra,dec = -90,-90
#     try:
#         header = fits.getheader(fits_file)
#         w = WCS(header)
#         ra,dec = w.wcs_pix2world(xcoord,ycoord,1)
#         keep_going = True
#     except Exception as e:
#         print('Function failed: get_radec_from_pixel_coords()')
#         print(e)
#         send_status_email('Function failed: mcd.get_radec_from_pixel_coords()','Function failed: mcd.get_radec_from_pixel_coords()')
#         keep_going = False
#     return ra,dec
# 
# 
# ##### FUNCTION DEFINITIONS -- DATABASE FUNCTIONS #####
# 
# def create_connection(db_file):
#     # create a database connection to the SQLite database specified by the db_file
#     # param db_file: database file
#     # return: Connection object or None
#     try:
#         conn = sqlite3.connect(db_file)
#         conn.text_factory = str
#         return conn
#     except Error as e:
#         print('Function failed: create_connection()')
#         print(e)
#     return None
# 
# def add_apostrophe_escape(ast_name):
#     # add escape sequence for apostrophes for use in SQL insertion command
#     # param ast_name: original asteroid name with unescaped apostrophe(s)
#     # return: ast_name_new: modified asteroid name with escaped apostrophe(s)
#     num_apostrophes = 0
#     name_length = len(ast_name)
#     for idx in range(0,name_length):
#         if ast_name[idx] == "'":
#             num_apostrophes += 1
#     ast_name_temp = [0 for idx in range(name_length+num_apostrophes)] # initialize new name string
#     idx_new = 0
#     for idx in range(0,name_length):
#         if ast_name[idx] != "'":
#             ast_name_temp[idx_new] = ast_name[idx]
#         else:
#             ast_name_temp[idx_new]   = "'"
#             ast_name_temp[idx_new+1] = "'"
#             idx_new += 1
#         idx_new += 1
#     ast_name_new = ''.join(ast_name_temp)
#     return ast_name_new
# 
# 
# ##### FUNCTION DEFINITIONS -- EMAILS #####
# 
# def send_status_email(subject,message_text):
#     port = 465  # For SSL
#     smtp_server = "smtp.gmail.com"
#     sender_email = "macadamia.asteroid.archive@gmail.com"  # Enter your address
#     receiver_email = "hhsieh@gmail.com"  # Enter receiver address
#     recipient = "Henry Hsieh <hhsieh@gmail.com>"
#     #password = input("Type your password and press enter: ")
#     password = "D0n'teatthecrabd1p!"
#     context = ssl.create_default_context()
#     message = '{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),message_text)
#     msg = 'To: {:s}\nSubject: {:s}\n\n{:s}'.format(recipient,subject,message)
#     #with smtplib.SMTP_SSL(smtp_server, port, context=context) as server_ssl:
#     #    server_ssl.login(sender_email,password)
#     #    server_ssl.sendmail(sender_email,receiver_email,msg)
#     return None
# 
# def send_status_email_old(subject,message):
#     try:
#         port = 465  # For SSL
#         smtp_server = "smtp.gmail.com"
#         sender_email = "macadamia.asteroid.archive@gmail.com"  # Enter your address
#         receiver_email = "hhsieh@gmail.com"  # Enter receiver address
#         recipient = "Henry Hsieh <hhsieh@gmail.com>"
#         #password = input("Type your password and press enter: ")
#         password = "D0n'teatthecrabd1p!"
#         context = ssl.create_default_context()
#         msg = 'To: {:s}\nSubject: {:s}\n\n{:s}'.format(recipient,subject,message)
#         with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
#             server.login(sender_email,password)
#             server.sendmail(sender_email,receiver_email,msg)
#     except Exception:
#         print('{:s} - Function failed: send_status_email()'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
#     return None
# 
# 
# ##### FUNCTION DEFINITIONS -- INITIALIZING FUNCTIONS #####
# 
# def validate_input_params(base_path,sqlite_file):
#     keep_going = True
#     if base_path[-1:] != '/': base_path = base_path + '/'
#     if not os.path.isdir(base_path):
#         print('Directory {:s} not found.'.format(base_path))
#         keep_going = False
#     if not os.path.exists(sqlite_file):
#         print('Database {:s} not found.'.format(sqlite_file))
#         keep_going = False
#     return base_path,keep_going
# 
# 
# def initialize_log_error_file(base_path,function_name):
#     # Create and open log file
#     dir_logfiles = base_path + 'log_files/'
#     if not os.path.isdir(dir_logfiles): os.mkdir(dir_logfiles)
#     path_logfile   = dir_logfiles + 'log_{:s}_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),function_name)
#     path_errorfile = dir_logfiles + 'errors_{:s}_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),function_name)
#     with open(path_logfile,'w') as log_file:
#         log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#     with open(path_errorfile,'w') as log_file:
#         log_file.write('Archival Asteroid Photometry Error Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#     print('Log file path:   {:s}'.format(path_logfile))
#     print('Error file path: {:s}'.format(path_errorfile))
#     return path_logfile,path_errorfile
# 
# 
# 
# def output_log_entry_nonstring(path_logfile,log_entry):
#     with open(path_logfile,'a') as log_file:
#         print('{:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(log_entry)
#         log_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(log_entry,file=log_file)
#     return None
# 
# def output_error_log_entry(path_logfile,path_errorfile,log_entry):
#     with open(path_logfile,'a') as log_file, open(path_errorfile,'a') as error_file:
#         print('{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
#         log_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
#         error_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
#     return None
# 
# def output_error_log_entry_nonstring(path_logfile,path_errorfile,log_entry):
#     with open(path_logfile,'a') as log_file, open(path_errorfile,'a') as error_file:
#         print('{:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(log_entry)
#         log_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(log_entry,file=log_file)
#         error_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(log_entry,file=error_file)
#     return None
# 
# def remove_error_log_if_empty(path_logfile,path_errorfile):
#     if os.path.exists(path_logfile) and os.path.exists(path_errorfile):
#         if len(open(path_errorfile).readlines()) < 2:
#             os.remove(path_errorfile)
#             output_log_entry(path_logfile,'Removing empty error log file.')
#     return None
# 
# def remove_output_file_if_empty(path_logfile,path_outputfile):
#     if os.path.exists(path_logfile) and os.path.exists(path_outputfile):
#         if len(open(path_outputfile).readlines()) < 2:
#             os.remove(path_outputfile)
#             output_log_entry(path_logfile,'Removing empty output file: {:s}'.format(path_outputfile))
#     return None
# 
# 
# ##### FUNCTION DEFINITIONS -- WCS FUNCTIONS #####
# 
# def get_pixel_coords_from_radec(fits_file,ra,dec):
#     px,py = [-1,-1],[-1,-1]
#     try:
#         header = fits.getheader(fits_file)
#         w = WCS(header)
#         px,py = w.wcs_world2pix(ra,dec,1)
#     except Exception as e:
#         print('Function failed: get_pixel_coords_from_radec()')
#         print(e)
#     return px,py
# 
# 
# def get_radec_from_pixel_coords(fits_file,xcoord,ycoord):
#     ra,dec = -90,-90
#     try:
#         header = fits.getheader(fits_file)
#         w = WCS(header)
#         ra,dec = w.wcs_pix2world(xcoord,ycoord,1)
#         keep_going = True
#     except Exception as e:
#         print('Function failed: get_radec_from_pixel_coords()')
#         print(e)
#         keep_going = False
#     return ra,dec
# 
# 
# def compute_angular_dist(ra1,dec1,ra2,dec2):
#     # Compute angular distance in degrees from two sets of RA,Dec coordinates, also in degrees
#     source_distance_deg = -1
#     try:
#         c1 = SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
#         c2 = SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
#         source_distance_deg = c1.separation(c2).degree
#     except Exception as e:
#         print('Function failed: compute_angular_dist()')
#         print(e)
#     return source_distance_deg
# 
# 
# def magavg(mag_array,magerr_array):
#     avgmag,avgmagerr = -999,-999
#     try:
#         num_mags = len(mag_array)
#         intensity_total = 0
#         error_total = 0
#         # initialize intensity array
#         intensity = [[0 for idx1 in range(2)] for idx2 in range(num_mags)]
#         for idx in range(0,num_mags):
#             intensity[idx][0] = 10**(0.4*(0.-mag_array[idx]))
#             intensity[idx][1] = ((magerr_array[idx]**2) * ((-0.921*(10**(-0.4*mag_array[idx])))**2.))**0.5
#             intensity_total   = intensity_total + intensity[idx][0] / (intensity[idx][1]**2)
#             error_total       = error_total + (intensity[idx][1]**2)**(-1)
#         avgintensity = intensity_total / error_total
#         avgintensityerr = error_total ** (-0.5)
#         avgmag = -2.5 * math.log10(avgintensity)
#         avgmagerr = 1.08574 * avgintensityerr/avgintensity
#     except Exception as e:
#         print('{:s} - Function failed: magavg()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
#         print(e)
#     return(avgmag,avgmagerr)
# 
# 
