import re
import sys
import time
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
import scipy
from scipy.stats import norm
from scipy import optimize
import urllib.request
from astropy.table import Table
import glob, os, bz2, subprocess
import os.path
from numpy import linspace
from decimal import *
import statistics
from astroquery.jplsbdb import SBDB
import macadamia_functions as mcd

##### UPDATED 5/6/21 #####


##### FUNCTION DEFINITIONS -- SDSS DATA HANDLING #####

def parse_sdssmoc_table(sdssmoc_filepath,sdssmoc_output_filepath,resolved_sdssmoc_name_log,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Reading in SDSS MOC data...')
    
            with open(sdssmoc_output_filepath,'w') as of:
                of.write('AstNum          MJD          RA         Dec   gmag gerr   rmag rerr   imag ierr   zmag zerr\n')
            with open(resolved_sdssmoc_name_log,'w') as of:
                of.write('Log of resolved SDSS MOC asteroid designations\n')
            
            num_lines = sum(1 for line in open(sdssmoc_filepath))
            print('>>> Scanning {:d} lines...'.format(num_lines))
    
            num_objs = 0
            print('>>> 0...',end='')
            with open(sdssmoc_filepath) as input_file:
                #for _ in range(2): #skip first two header lines
                #    next(input_file)
                for line in input_file:
                    num_objs += 1
                    if num_objs % 10000  == 0: print('{:d}...'.format(num_objs),end='')
                    if num_objs % 100000 == 0: print('\n>>> ',end='')
                    mjd       = float(line[46:57])
                    ra        = float(line[58:68])
                    dec       = float(line[70:79])
                    gmag      = float(line[174:179])
                    gerr      = float(line[180:184])
                    rmag      = float(line[185:190])
                    rerr      = float(line[191:195])
                    imag      = float(line[196:201])
                    ierr      = float(line[202:206])
                    zmag      = float(line[207:212])
                    zerr      = float(line[213:217])
                    flag_id   = int(line[242:243])
                    ast_num   = int(line[245:251])
                    #ast_desig = int(line[245:251])
                    prov_desig_temp = line[252:264].strip()
                    if ast_num == 0 and prov_desig_temp != '-':
                        #with open(resolved_sdssmoc_name_log,'a') as of:
                        #     of.write('{:s} - {:s} ==> {:d}.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),prov_desig,ast_num))
                        prov_desig = ''
                        for idx in range(0,len(prov_desig_temp)):
                            if prov_desig_temp[idx] != '_': prov_desig = prov_desig + prov_desig_temp[idx]
                            else: prov_desig = prov_desig + ' '
                        sbdb_data = SBDB.query(prov_desig)
                        if 'object' in sbdb_data:
                            if sbdb_data['object']['des'].isnumeric():
                                ast_num = int(sbdb_data['object']['des'])
                                with open(resolved_sdssmoc_name_log,'a') as of:
                                    of.write('{:s} - {:s} ==> {:d}.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),prov_desig,ast_num))
                            else:
                                with open(resolved_sdssmoc_name_log,'a') as of:
                                    of.write('{:s} - No numerical designation available for {:s}.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),prov_desig))
                        else:
                            with open(resolved_sdssmoc_name_log,'a') as of:
                                of.write('{:s} - JPL look-up of {:s} failed.\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),prov_desig))
                    if flag_id == 1 and ast_num != 0:
                        with open(sdssmoc_output_filepath,'a') as of:
                            of.write('{:>6d}  {:11.5f}  {:10.6f}  {:10.6f}  {:5.2f} {:4.2f}  {:5.2f} {:4.2f}  {:5.2f} {:4.2f}  {:5.2f} {:4.2f}\n'.format(ast_num,mjd,ra,dec,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr))
                        #sdssmoc_table.add_row((ast_num,mjd,ra,dec,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr))

            mcd.output_log_entry(path_logfile,'\n Parsing SDSS MOC data complete.')
        
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: parse_sdssmoc_table()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:    
        mcd.output_log_entry(path_logfile,'Function skipped: parse_sdssmoc_table()')
    return keep_going


def create_sdssmoc_table(sdssmoc_output_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Creating SDSS MOC table...')
    
            ast_num = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(0)) 
            mjd     = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(1)) 
            ra      = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(2)) 
            dec     = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(3)) 
            gmag    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(4)) 
            gerr    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(5)) 
            rmag    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(6)) 
            rerr    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(7)) 
            imag    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(8)) 
            ierr    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(9)) 
            zmag    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(10)) 
            zerr    = np.genfromtxt(sdssmoc_output_filepath,skip_header=1,usecols=(11)) 
    
            mcd.output_log_entry(path_logfile,'>>> SDSS MOC table created with {:d} detections of numbered asteroids...'.format(len(ast_num)))
    
            sdssmoc_table = Table([ast_num,mjd,ra,dec,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr],
                                    names=('ast_num','mjd','ra','dec','gmag','gerr','rmag','rerr','imag','ierr','zmag','zerr'), \
                                    dtype=('i8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

            sdssmoc_table.sort('ast_num')  # sort SDSS MOC table in increasing asteroid numerical designation
    
            mcd.output_log_entry(path_logfile,'Creation of SDSS MOC table complete.')
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_sdssmoc_table()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:    
        mcd.output_log_entry(path_logfile,'Function skipped: create_sdssmoc_table()')
    return sdssmoc_table,keep_going


def create_sdss_macadamia_comparison_files(sdssmoc_table,macadamia_output_filepath,sdss_vs_macadamia_filepath,keep_going,path_logfile,path_errorfile):
    outlier_threshold       = 0.3  # maximum allowable abs(mag_macadamia-mag_sdss) for outlier removal
    highprecision_threshold = 0.1  # maximum mag err for either SDSS or MACADAMIA photometry for high-precision data selection
    if keep_going:
        mcd.output_log_entry(path_logfile,'Creating SDSS/MACADAMIA comparison files...')
        #mcd.decompress_file_gzip(macadamia_output_filepath)
        macadamia_output_filepath       = macadamia_output_filepath[:-3]
        sdss_vs_macadamia_filepath_all  = sdss_vs_macadamia_filepath[:-4] + '_all.dat'
        sdss_vs_macadamia_filepath_g    = sdss_vs_macadamia_filepath[:-4] + '_g.dat'
        sdss_vs_macadamia_filepath_r    = sdss_vs_macadamia_filepath[:-4] + '_r.dat'
        sdss_vs_macadamia_filepath_i    = sdss_vs_macadamia_filepath[:-4] + '_i.dat'
        sdss_vs_macadamia_filepath_z    = sdss_vs_macadamia_filepath[:-4] + '_z.dat'
        sdss_vs_macadamia_filepath_g_or = sdss_vs_macadamia_filepath[:-4] + '_g_or.dat'
        sdss_vs_macadamia_filepath_r_or = sdss_vs_macadamia_filepath[:-4] + '_r_or.dat'
        sdss_vs_macadamia_filepath_i_or = sdss_vs_macadamia_filepath[:-4] + '_i_or.dat'
        sdss_vs_macadamia_filepath_z_or = sdss_vs_macadamia_filepath[:-4] + '_z_or.dat'
        sdss_vs_macadamia_filepath_g_hp = sdss_vs_macadamia_filepath[:-4] + '_g_hp.dat'
        sdss_vs_macadamia_filepath_r_hp = sdss_vs_macadamia_filepath[:-4] + '_r_hp.dat'
        sdss_vs_macadamia_filepath_i_hp = sdss_vs_macadamia_filepath[:-4] + '_i_hp.dat'
        sdss_vs_macadamia_filepath_z_hp = sdss_vs_macadamia_filepath[:-4] + '_z_hp.dat'
                        
        mcd.output_log_entry(path_logfile,'Reading in MACADAMIA data...')
        mcd.output_log_entry(path_logfile,'Reading in detection IDs...')
        detection_id_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(0))
        mcd.output_log_entry(path_logfile,'Reading in numerical designations...')
        desig_num_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(1))
        mcd.output_log_entry(path_logfile,'Reading in exposure start times...')
        exp_start_jd_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(10))
        mcd.output_log_entry(path_logfile,'Reading in distances to predicted positions...')
        dist_pred_postn_mac = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(31))
        mcd.output_log_entry(path_logfile,'Reading in measured magnitudes and uncertainties...')
        mag_mac             = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(34))
        magerr_mac          = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(35))
        mcd.output_log_entry(path_logfile,'Reading in SNRs...')
        snr_mac             = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(37))
        mcd.output_log_entry(path_logfile,'Reading in trail lengths and FWHMs...')
        trail_length_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(38))
        trail_fwhm_mac      = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(39))
        trail_len_pred_mac  = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(71))
        mcd.output_log_entry(path_logfile,'Reading in distances and magnitudes to closest sources...')
        dist_src1_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(44))
        mag_src1_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(45))
        dist_src2_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(50))
        mag_src2_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(51))
        dist_src3_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(56))
        mag_src3_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(57))
        mcd.output_log_entry(path_logfile,'Reading in global and local source densities...')
        src_density_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(21))
        local_srcdn_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(59))
        mcd.output_log_entry(path_logfile,'Reading in distances to edges...')
        dist_edge_l_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(60))
        dist_edge_r_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(61))
        dist_edge_b_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(62))
        dist_edge_t_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(63))
        mcd.output_log_entry(path_logfile,'Reading in zero point parameters...')
        zpoint_mac          = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(15))
        zpt_nstars_mac      = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(17))
        mcd.output_log_entry(path_logfile,'Reading in limiting magnitudes...')
        limit_mag_ps_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(18))
        limit_mag_sb_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(19))
        mcd.output_log_entry(path_logfile,'Reading in seeing measurements...')
        psf_width_mean_mac  = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(20))
        psf_width_src_mac   = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(36))
        num_dets_mac        = len(detection_id_mac)

        mcd.output_log_entry(path_logfile,'Reading in filter names...')
        filter_name_mac  = [0 for idx in range(num_dets_mac)]
        idx_mac = 0
        with open(macadamia_output_filepath) as input_file:
            for _ in range(1): #skip first header line
                next(input_file)
            for line in input_file:
                filter_name_mac[idx_mac]  = line[125:135].strip()
                idx_mac += 1

        with open(sdss_vs_macadamia_filepath_all,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')

        outliers_g,outliers_r,outliers_i,outliers_z = 0,0,0,0
        higherr_g,higherr_r,higherr_i,higherr_z = 0,0,0,0
        for idx_mac in range(0,num_dets_mac):
            idx_sdss = 0
            sdss_entry_found = False
            
            # find distance to nearest chip edge
            dist_edge_mac = dist_edge_l_mac[idx_mac]
            if dist_edge_r_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_r_mac[idx_mac]
            if dist_edge_t_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_t_mac[idx_mac]
            if dist_edge_b_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_b_mac[idx_mac]
                
            mcd.output_log_entry(path_logfile,'Testing MACADAMIA detection {:d}: asteroid {:d},{:f},mag={:f}'.format(idx_mac,desig_num_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac]))
            if mag_mac[idx_mac] != 99.99:
                while idx_sdss < len(sdssmoc_table) and not sdss_entry_found:
                    #mcd.output_log_entry(path_logfile,'SDSS detection {:d}: asteroid {:d},{:f} :: MACADAMIA detection {:d}: asteroid {:d},{:f}'.format(idx_sdss,sdssmoc_table['ast_num'][idx_sdss],sdssmoc_table['mjd'][idx_sdss],idx_mac,desig_num_mac[idx_mac],exp_start_jd_mac[idx_mac]))
                    if desig_num_mac[idx_mac] == sdssmoc_table['ast_num'][idx_sdss]:
                        if abs((exp_start_jd_mac[idx_mac]-2400000.5)-sdssmoc_table['mjd'][idx_sdss]) < 0.005:
                            sdss_entry_found = True
                            mcd.output_log_entry(path_logfile,'SDSS entry {:d} matched to detection_id {:d}.'.format(idx_sdss,int(detection_id_mac[idx_mac])))
                            if filter_name_mac[idx_mac] == 'g' and sdssmoc_table['gmag'][idx_sdss] != 99.99:
                                mag_sdss          = sdssmoc_table['gmag'][idx_sdss]
                                merr_sdss         = sdssmoc_table['gerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'r' and sdssmoc_table['rmag'][idx_sdss] != 99.99:
                                mag_sdss          = sdssmoc_table['rmag'][idx_sdss]
                                merr_sdss         = sdssmoc_table['rerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'i' and sdssmoc_table['imag'][idx_sdss] != 99.99:
                                mag_sdss          = sdssmoc_table['imag'][idx_sdss]
                                merr_sdss         = sdssmoc_table['ierr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'z' and sdssmoc_table['zmag'][idx_sdss] != 99.99:
                                mag_sdss          = sdssmoc_table['zmag'][idx_sdss]
                                merr_sdss         = sdssmoc_table['zerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            mcd.output_log_entry(path_logfile,'> SDSS entry {:d} removed; {:d} SDSS entries remaining.'.format(idx_sdss,len(sdssmoc_table)))
                            sdssmoc_table.remove_row(idx_sdss)
                    idx_sdss += 1
            
            if sdss_entry_found:
                # write detection data to sdss vs. macadamia output file
                if snr_mac[idx_mac] != -999:
                    with open(sdss_vs_macadamia_filepath_all,'a') as of:
                        of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                    if filter_name_mac[idx_mac] == 'g':                            
                        with open(sdss_vs_macadamia_filepath_g,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                            with open(sdss_vs_macadamia_filepath_g_or,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            outliers_g += 1
                        if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                            with open(sdss_vs_macadamia_filepath_g_hp,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            higherr_g += 1
                    if filter_name_mac[idx_mac] == 'r':
                        with open(sdss_vs_macadamia_filepath_r,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                            with open(sdss_vs_macadamia_filepath_r_or,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            outliers_r += 1
                        if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                            with open(sdss_vs_macadamia_filepath_r_hp,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            higherr_r += 1
                    if filter_name_mac[idx_mac] == 'i':
                        with open(sdss_vs_macadamia_filepath_i,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                            with open(sdss_vs_macadamia_filepath_i_or,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            outliers_i += 1
                        if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                            with open(sdss_vs_macadamia_filepath_i_hp,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            higherr_i += 1
                    if filter_name_mac[idx_mac] == 'z':
                        with open(sdss_vs_macadamia_filepath_z,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                            with open(sdss_vs_macadamia_filepath_z_or,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            outliers_z += 1
                        if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                            with open(sdss_vs_macadamia_filepath_z_hp,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        else:
                            higherr_z += 1
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'SDSS match for detection_id {:d} not found.'.format(int(detection_id_mac[idx_mac])))
                if snr_mac[idx_mac] != -999:
                    with open(sdss_vs_macadamia_filepath,'a') as of:
                        of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],0,-9.99,-9.99,-9.99,-9.99,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from g'-band photometry.".format(outliers_g,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from r'-band photometry.".format(outliers_r,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from i'-band photometry.".format(outliers_i,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from z'-band photometry.".format(outliers_z,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from g'-band photometry.".format(higherr_g,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from r'-band photometry.".format(higherr_r,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from i'-band photometry.".format(higherr_i,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from z'-band photometry.".format(higherr_z,highprecision_threshold))
        #mcd.compress_file_gzip(macadamia_output_filepath)
        mcd.output_log_entry(path_logfile,'Creation of SDSS/MACADAMIA comparison files complete.')
    else:    
        mcd.output_log_entry(path_logfile,'Function skipped: create_sdss_macadamia_comparison_files()')
    
    return keep_going


def create_sdss_macadamia_comparison_files2(sdssmoc_table,macadamia_output_filepath,sdss_vs_macadamia_filepath,keep_going,path_logfile,path_errorfile):
    outlier_threshold       = 0.3  # maximum allowable abs(mag_macadamia-mag_sdss) for outlier removal
    highprecision_threshold = 0.1  # maximum mag err for either SDSS or MACADAMIA photometry for high-precision data selection
    if keep_going:
        mcd.output_log_entry(path_logfile,'Creating SDSS/MACADAMIA comparison files...')
        #mcd.decompress_file_gzip(macadamia_output_filepath)
        macadamia_output_filepath           = macadamia_output_filepath[:-3]
        sdss_vs_macadamia_filepath_all      = sdss_vs_macadamia_filepath[:-4] + '_all.dat'
        sdss_vs_macadamia_filepath_notfound = sdss_vs_macadamia_filepath[:-4] + '_notfound.dat'
        sdss_vs_macadamia_filepath_g        = sdss_vs_macadamia_filepath[:-4] + '_g.dat'
        sdss_vs_macadamia_filepath_r        = sdss_vs_macadamia_filepath[:-4] + '_r.dat'
        sdss_vs_macadamia_filepath_i        = sdss_vs_macadamia_filepath[:-4] + '_i.dat'
        sdss_vs_macadamia_filepath_z        = sdss_vs_macadamia_filepath[:-4] + '_z.dat'
        sdss_vs_macadamia_filepath_g_or     = sdss_vs_macadamia_filepath[:-4] + '_g_or.dat'
        sdss_vs_macadamia_filepath_r_or     = sdss_vs_macadamia_filepath[:-4] + '_r_or.dat'
        sdss_vs_macadamia_filepath_i_or     = sdss_vs_macadamia_filepath[:-4] + '_i_or.dat'
        sdss_vs_macadamia_filepath_z_or     = sdss_vs_macadamia_filepath[:-4] + '_z_or.dat'
        sdss_vs_macadamia_filepath_g_hp     = sdss_vs_macadamia_filepath[:-4] + '_g_hp.dat'
        sdss_vs_macadamia_filepath_r_hp     = sdss_vs_macadamia_filepath[:-4] + '_r_hp.dat'
        sdss_vs_macadamia_filepath_i_hp     = sdss_vs_macadamia_filepath[:-4] + '_i_hp.dat'
        sdss_vs_macadamia_filepath_z_hp     = sdss_vs_macadamia_filepath[:-4] + '_z_hp.dat'
                        
        mcd.output_log_entry(path_logfile,'Reading in MACADAMIA data...')
        mcd.output_log_entry(path_logfile,'Reading in detection IDs...')
        detection_id_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(0))
        mcd.output_log_entry(path_logfile,'Reading in numerical designations...')
        desig_num_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(1))
        mcd.output_log_entry(path_logfile,'Reading in exposure start times...')
        exp_start_jd_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(10))
        mcd.output_log_entry(path_logfile,'Reading in distances to predicted positions...')
        dist_pred_postn_mac = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(31))
        mcd.output_log_entry(path_logfile,'Reading in measured magnitudes and uncertainties...')
        mag_mac             = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(34))
        magerr_mac          = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(35))
        mcd.output_log_entry(path_logfile,'Reading in SNRs...')
        snr_mac             = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(37))
        mcd.output_log_entry(path_logfile,'Reading in trail lengths and FWHMs...')
        trail_length_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(38))
        trail_fwhm_mac      = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(39))
        trail_len_pred_mac  = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(71))
        mcd.output_log_entry(path_logfile,'Reading in distances and magnitudes to closest sources...')
        dist_src1_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(44))
        mag_src1_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(45))
        dist_src2_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(50))
        mag_src2_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(51))
        dist_src3_mac       = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(56))
        mag_src3_mac        = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(57))
        mcd.output_log_entry(path_logfile,'Reading in global and local source densities...')
        src_density_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(21))
        local_srcdn_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(59))
        mcd.output_log_entry(path_logfile,'Reading in distances to edges...')
        dist_edge_l_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(60))
        dist_edge_r_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(61))
        dist_edge_b_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(62))
        dist_edge_t_mac     = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(63))
        mcd.output_log_entry(path_logfile,'Reading in zero point parameters...')
        zpoint_mac          = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(15))
        zpt_nstars_mac      = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(17))
        mcd.output_log_entry(path_logfile,'Reading in limiting magnitudes...')
        limit_mag_ps_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(18))
        limit_mag_sb_mac    = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(19))
        mcd.output_log_entry(path_logfile,'Reading in seeing measurements...')
        psf_width_mean_mac  = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(20))
        psf_width_src_mac   = np.genfromtxt(macadamia_output_filepath,skip_header=1,usecols=(36))
        num_dets_mac        = len(detection_id_mac)

        mcd.output_log_entry(path_logfile,'Reading in filter names...')
        filter_name_mac  = [0 for idx in range(num_dets_mac)]
        idx_mac = 0
        with open(macadamia_output_filepath) as input_file:
            for _ in range(1): #skip first header line
                next(input_file)
            for line in input_file:
                filter_name_mac[idx_mac]  = line[125:135].strip()
                idx_mac += 1

        with open(sdss_vs_macadamia_filepath_notfound,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_all,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z_or,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_g_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_r_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_i_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')
        with open(sdss_vs_macadamia_filepath_z_hp,'w') as of:
            of.write('detection_id   desig_num   filt       mjd_start     mag   merr      mjd_sdss   mag_sdss  merr_sdss   mag_diff  mag_diff_err       snr   src_density   local_srcdn   dist_pred_postn   trail_predln   trail_length   trail_fwhm   dist_src1  mag_src1   dist_src2  mag_src2   dist_src3  mag_src3   dist_edge   zpoint   zpt_nstars   lim_mag_ps   lim_mag_sb   mean_seeing   src_seeing\n')

        outliers_g,outliers_r,outliers_i,outliers_z = 0,0,0,0
        higherr_g,higherr_r,higherr_i,higherr_z = 0,0,0,0
        
        for idx_sdss in range(0,len(sdssmoc_table)):
            
            mcd.output_log_entry(path_logfile,'Matching asteroid {:d} at {:.5f} (SDSS entry {:d})...'.format(int(sdssmoc_table['ast_num'][idx_sdss]),sdssmoc_table['mjd'][idx_sdss],idx_sdss))
            sdss_filters = ['g','r','i','z']
            for sdss_filter in sdss_filters:
                #print(' >>> {:s}-band...'.format(sdss_filter))
                #log_file.write(' >>> {:s}-band...\n'.format(sdss_filter))
                if sdss_filter == 'g':
                    mag_sdss  = sdssmoc_table['gmag'][idx_sdss]
                    merr_sdss = sdssmoc_table['gerr'][idx_sdss]
                if sdss_filter == 'r':
                    mag_sdss  = sdssmoc_table['rmag'][idx_sdss]
                    merr_sdss = sdssmoc_table['rerr'][idx_sdss]
                if sdss_filter == 'i':
                    mag_sdss  = sdssmoc_table['imag'][idx_sdss]
                    merr_sdss = sdssmoc_table['ierr'][idx_sdss]
                if sdss_filter == 'z':
                    mag_sdss  = sdssmoc_table['zmag'][idx_sdss]
                    merr_sdss = sdssmoc_table['zerr'][idx_sdss]
                idx_mac = 0
                mac_entry_found = False
                while idx_mac < len(detection_id_mac) and not mac_entry_found:
                    #mcd.output_log_entry(path_logfile,'Testing MACADAMIA detection {:d}: asteroid {:d} at {:f}'.format(idx_mac,int(desig_num_mac[idx_mac]),exp_start_jd_mac[idx_mac]))
                    if desig_num_mac[idx_mac] == sdssmoc_table['ast_num'][idx_sdss]:
                        if abs((exp_start_jd_mac[idx_mac]-2400000.5)-sdssmoc_table['mjd'][idx_sdss]) < 0.005 and filter_name_mac[idx_mac] == sdss_filter:
                            mac_entry_found = True
                            mcd.output_log_entry(path_logfile,'>> SDSS entry {:d} ({:s}-band) matched to MACADAMIA detection {:d}.'.format(idx_sdss,sdss_filter,int(detection_id_mac[idx_mac])))
                            if filter_name_mac[idx_mac] == 'g' and sdssmoc_table['gmag'][idx_sdss] != 99.99:
                                #mag_sdss          = sdssmoc_table['gmag'][idx_sdss]
                                #merr_sdss         = sdssmoc_table['gerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'r' and sdssmoc_table['rmag'][idx_sdss] != 99.99:
                                #mag_sdss          = sdssmoc_table['rmag'][idx_sdss]
                                #merr_sdss         = sdssmoc_table['rerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'i' and sdssmoc_table['imag'][idx_sdss] != 99.99:
                                #mag_sdss          = sdssmoc_table['imag'][idx_sdss]
                                #merr_sdss         = sdssmoc_table['ierr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                            if filter_name_mac[idx_mac] == 'z' and sdssmoc_table['zmag'][idx_sdss] != 99.99:
                                #mag_sdss          = sdssmoc_table['zmag'][idx_sdss]
                                #merr_sdss         = sdssmoc_table['zerr'][idx_sdss]
                                mag_diff_sdss_mac = mag_mac[idx_mac] - mag_sdss
                                mag_diff_err_sdss_mac = (magerr_mac[idx_mac]**2 + merr_sdss**2)**0.5
                    idx_mac += 1
                
                if mac_entry_found:
                    idx_mac -= 1
                    # find distance to nearest chip edge
                    dist_edge_mac = dist_edge_l_mac[idx_mac]
                    if dist_edge_r_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_r_mac[idx_mac]
                    if dist_edge_t_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_t_mac[idx_mac]
                    if dist_edge_b_mac[idx_mac] < dist_edge_mac: dist_edge_mac = dist_edge_b_mac[idx_mac]

                    # write detection data to sdss vs. macadamia output file
                    if snr_mac[idx_mac] != -999:
                        with open(sdss_vs_macadamia_filepath_all,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                        if filter_name_mac[idx_mac] == 'g':                            
                            with open(sdss_vs_macadamia_filepath_g,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                                with open(sdss_vs_macadamia_filepath_g_or,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                outliers_g += 1
                            if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                                with open(sdss_vs_macadamia_filepath_g_hp,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                higherr_g += 1
                        if filter_name_mac[idx_mac] == 'r':
                            with open(sdss_vs_macadamia_filepath_r,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                                with open(sdss_vs_macadamia_filepath_r_or,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                outliers_r += 1
                            if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                                with open(sdss_vs_macadamia_filepath_r_hp,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                higherr_r += 1
                        if filter_name_mac[idx_mac] == 'i':
                            with open(sdss_vs_macadamia_filepath_i,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                                with open(sdss_vs_macadamia_filepath_i_or,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                outliers_i += 1
                            if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                                with open(sdss_vs_macadamia_filepath_i_hp,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                higherr_i += 1
                        if filter_name_mac[idx_mac] == 'z':
                            with open(sdss_vs_macadamia_filepath_z,'a') as of:
                                of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            if abs(mag_mac[idx_mac] - mag_sdss) < outlier_threshold:
                                with open(sdss_vs_macadamia_filepath_z_or,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                outliers_z += 1
                            if magerr_mac[idx_mac] < highprecision_threshold and merr_sdss < highprecision_threshold:
                                with open(sdss_vs_macadamia_filepath_z_hp,'a') as of:
                                    of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(int(detection_id_mac[idx_mac]),int(desig_num_mac[idx_mac]),filter_name_mac[idx_mac],exp_start_jd_mac[idx_mac],mag_mac[idx_mac],magerr_mac[idx_mac],sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,mag_diff_sdss_mac,mag_diff_err_sdss_mac,snr_mac[idx_mac],src_density_mac[idx_mac],local_srcdn_mac[idx_mac],dist_pred_postn_mac[idx_mac],trail_len_pred_mac[idx_mac],trail_length_mac[idx_mac],trail_fwhm_mac[idx_mac],dist_src1_mac[idx_mac],mag_src1_mac[idx_mac],dist_src2_mac[idx_mac],mag_src2_mac[idx_mac],dist_src3_mac[idx_mac],mag_src3_mac[idx_mac],dist_edge_mac,zpoint_mac[idx_mac],int(zpt_nstars_mac[idx_mac]),limit_mag_ps_mac[idx_mac],limit_mag_sb_mac[idx_mac],psf_width_mean_mac[idx_mac],psf_width_src_mac[idx_mac]))
                            else:
                                higherr_z += 1
                else:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'>> MACADAMIA match for SDSS entry {:d} ({:s}-band) not found.'.format(idx_sdss,sdss_filter))
                    if mag_sdss != 99.99:
                        with open(sdss_vs_macadamia_filepath_notfound,'a') as of:
                            of.write('{:>12d}   {:>9d}      {:1s}   {:11.5f}   {:5.2f}  {:5.2f}   {:11.5f}      {:5.2f}      {:5.2f}      {:5.2f}         {:5.2f}   {:7.1f}   {:11.3f}   {:11.3f}   {:15.3f}   {:12.3f}   {:12.3f}   {:10.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.3f}  {:8.3f}   {:9.1f}   {:6.3f}   {:>10d}   {:10.3f}   {:10.3f}   {:11.3f}   {:10.3f}\n'.format(0,int(sdssmoc_table['ast_num'][idx_sdss]),sdss_filter,0,0,0,sdssmoc_table['mjd'][idx_sdss],mag_sdss,merr_sdss,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from g'-band photometry.".format(outliers_g,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from r'-band photometry.".format(outliers_r,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from i'-band photometry.".format(outliers_i,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} outliers (dev > {:1f} mag) removed from z'-band photometry.".format(outliers_z,outlier_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from g'-band photometry.".format(higherr_g,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from r'-band photometry.".format(higherr_g,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from i'-band photometry.".format(higherr_g,highprecision_threshold))
        mcd.output_log_entry(path_logfile,"{:d} low-precision data points (magerr_mcdm or magerr_sdss > {:1f} mag) removed from z'-band photometry.".format(higherr_g,highprecision_threshold))
        #mcd.compress_file_gzip(macadamia_output_filepath)
        mcd.output_log_entry(path_logfile,'Creation of SDSS/MACADAMIA comparison files complete.')
    else:    
        mcd.output_log_entry(path_logfile,'Function skipped: create_sdss_macadamia_comparison_files()')
    
    return keep_going

def create_sdss_macadamia_comparison_plots(sdss_vs_macadamia_filepath,dataset,keep_going,path_logfile,path_errorfile):
    if dataset == 'outliers_removed': filename_suffix = '_or'
    elif dataset == 'high_precision': filename_suffix = '_hp'
    else: filename_suffix = ''
    if keep_going:
        try:
            sdss_vs_macadamia_dev_vs_snr_all_filepath  = sdss_vs_macadamia_filepath[:-4] + '_dev_vs_snr_all.pdf'
            sdss_vs_macadamia_dev_vs_snr_griz_filepath = sdss_vs_macadamia_filepath[:-4] + '_dev_vs_snr_griz.pdf'
            sdss_vs_macadamia_dev_vs_mag_all_filepath  = sdss_vs_macadamia_filepath[:-4] + '_dev_vs_mag_all.pdf'
            sdss_vs_macadamia_dev_vs_mag_griz_filepath = sdss_vs_macadamia_filepath[:-4] + '_dev_vs_mag_griz.pdf'
    
            merr_g          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(5)) 
            merr_r          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(5)) 
            merr_i          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(5)) 
            merr_z          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(5)) 
            sdss_mags_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(7)) 
            sdss_mags_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(7)) 
            sdss_mags_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(7)) 
            sdss_mags_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(7)) 
            sdss_merr_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(8)) 
            sdss_merr_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(8)) 
            sdss_merr_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(8)) 
            sdss_merr_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(8)) 
            mag_diffs_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(9)) 
            mag_diffs_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(9)) 
            mag_diffs_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(9)) 
            mag_diffs_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(9)) 
            mag_diff_errs_g = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(10)) 
            mag_diff_errs_r = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(10)) 
            mag_diff_errs_i = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(10)) 
            mag_diff_errs_z = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(10)) 
            snrs_g          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(11)) 
            snrs_r          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(11)) 
            snrs_i          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(11)) 
            snrs_z          = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(11)) 
            src_density_g   = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(12)) 
            src_density_r   = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(12)) 
            src_density_i   = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(12)) 
            src_density_z   = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(12)) 
            dist_predpstn_g = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(13)) 
            dist_predpstn_r = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(13)) 
            dist_predpstn_i = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(13)) 
            dist_predpstn_z = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(13)) 
            trail_length_g  = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(14)) 
            trail_length_r  = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(14)) 
            trail_length_i  = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(14)) 
            trail_length_z  = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(14)) 
            trail_fwhm_g    = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(15)) 
            trail_fwhm_r    = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(15)) 
            trail_fwhm_i    = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(15)) 
            trail_fwhm_z    = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(15)) 
            dist_src1_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(16)) 
            dist_src1_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(16)) 
            dist_src1_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(16)) 
            dist_src1_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(16)) 
            mag_src1_g      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(17)) 
            mag_src1_r      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(17)) 
            mag_src1_i      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(17)) 
            mag_src1_z      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(17)) 
            dist_src2_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(18)) 
            dist_src2_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(18)) 
            dist_src2_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(18)) 
            dist_src2_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(18)) 
            mag_src2_g      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(19)) 
            mag_src2_r      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(19)) 
            mag_src2_i      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(19)) 
            mag_src2_z      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(19)) 
            dist_src3_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(20)) 
            dist_src3_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(20)) 
            dist_src3_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(20)) 
            dist_src3_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(20)) 
            mag_src3_g      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(21)) 
            mag_src3_r      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(21)) 
            mag_src3_i      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(21)) 
            mag_src3_z      = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(21)) 
            dist_edge_g     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_g'+filename_suffix+'.dat',skip_header=1,usecols=(22)) 
            dist_edge_r     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_r'+filename_suffix+'.dat',skip_header=1,usecols=(22)) 
            dist_edge_i     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_i'+filename_suffix+'.dat',skip_header=1,usecols=(22)) 
            dist_edge_z     = np.genfromtxt(sdss_vs_macadamia_filepath[:-4]+'_z'+filename_suffix+'.dat',skip_header=1,usecols=(22)) 

            #keep_going,log_file = create_histogram_all(mag_diffs_all,sdss_vs_macadamia_filepath[:-4]+'_histogram_all.pdf',keep_going,log_file)
            #keep_going,log_file = create_dev_vs_mag_all(mag_diffs_all,mags_all,sdss_vs_macadamia_filepath[:-4]+'_dev_vs_mag_all.pdf',keep_going,log_file)
            #keep_going,log_file = create_dev_vs_snr_all(mag_diffs_all,snrs_all,sdss_vs_macadamia_filepath[:-4]+'_dev_vs_snr_all.pdf',keep_going,log_file)

            keep_going = create_err_histograms_griz(merr_g,merr_r,merr_i,merr_z,sdss_merr_g,sdss_merr_r,sdss_merr_i,sdss_merr_z,sdss_vs_macadamia_filepath[:-4]+'_err_histogram_griz'+filename_suffix+'.pdf',keep_going,path_logfile,path_errorfile)
            keep_going = create_diff_histograms_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,sdss_vs_macadamia_filepath[:-4]+'_histogram_griz'+filename_suffix+'.pdf',keep_going,path_logfile,path_errorfile)
            keep_going = create_dev_vs_mag_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,mags_g,mags_r,mags_i,mags_z,sdss_vs_macadamia_filepath[:-4]+'_dev_vs_mag_griz'+filename_suffix+'.pdf',keep_going,path_logfile,path_errorfile)
            keep_going = create_dev_vs_snr_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,snrs_g,snrs_r,snrs_i,snrs_z,sdss_vs_macadamia_filepath[:-4]+'_dev_vs_snr_griz'+filename_suffix+'.pdf',keep_going,path_logfile,path_errorfile)

        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_sdss_macadamia_comparison_plots()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:    
        mcd.output_log_entry(path_logfile,'Function skipped: create_sdss_macadamia_comparison_plots()')
    
    return keep_going

def gaussian(x,amplitude,mean,stddev):
    return amplitude * np.exp(-((x-mean)**2)/(2*stddev**2))

def create_err_histograms_griz(merr_g,merr_r,merr_i,merr_z,sdss_merr_g,sdss_merr_r,sdss_merr_i,sdss_merr_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1a = [0.12,0.63,0.35,0.20]
            rect_plot_1b = [0.12,0.40,0.35,0.20]
            rect_plot_2a = [0.12,0.32,0.35,0.20]
            rect_plot_2b = [0.12,0.12,0.35,0.20]
            rect_plot_3a = [0.55,0.63,0.35,0.20]
            rect_plot_3b = [0.55,0.40,0.35,0.20]
            rect_plot_4a = [0.55,0.32,0.35,0.20]
            rect_plot_4b = [0.55,0.12,0.35,0.20]
            plt.figure(1,figsize=(12,12))
            axPlot_g1 = plt.axes(rect_plot_1a)
            axPlot_g2 = plt.axes(rect_plot_1b)
            axPlot_r1 = plt.axes(rect_plot_2a)
            axPlot_r2 = plt.axes(rect_plot_2b)
            axPlot_i1 = plt.axes(rect_plot_3a)
            axPlot_i2 = plt.axes(rect_plot_3b)
            axPlot_z1 = plt.axes(rect_plot_4a)
            axPlot_z2 = plt.axes(rect_plot_4b)
            
            n_mag_err_bins = 20
            #mag_diff_min = np.amin(mag_diffs)
            #mag_diff_max = np.amax(mag_diffs)
            #n_mag_diff_bins = 401
            mag_err_min = 0
            mag_err_max = 1
            mag_err_bins = np.arange(mag_err_min,mag_err_max,(mag_err_max-mag_err_min)/n_mag_err_bins)
            #mu_g,sigma_g = norm.fit(mag_diffs_g)
            #mu_r,sigma_r = norm.fit(mag_diffs_r)
            #mu_i,sigma_i = norm.fit(mag_diffs_i)
            #mu_z,sigma_z = norm.fit(mag_diffs_z)
            x_g1,bins_g1,p_g1 = axPlot_g1.hist(merr_g,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_r1,bins_r1,p_r1 = axPlot_r1.hist(merr_r,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_i1,bins_i1,p_i1 = axPlot_i1.hist(merr_i,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_z1,bins_z1,p_z1 = axPlot_z1.hist(merr_z,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_g2,bins_g2,p_g2 = axPlot_g2.hist(sdss_merr_g,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_r2,bins_r2,p_r2 = axPlot_r2.hist(sdss_merr_r,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_i2,bins_i2,p_i2 = axPlot_i2.hist(sdss_merr_i,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_z2,bins_z2,p_z2 = axPlot_z2.hist(sdss_merr_z,bins=mag_err_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            
            major_xticks = np.arange(-2,2,0.2)
            minor_xticks = np.arange(-2,2,0.1)
            #major_xticks = np.arange(-12,12,2)
            #minor_xticks = np.arange(-12,12,0.5)
            major_yticks = np.arange(-2,2,0.1)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g1.tick_params(axis='both',which='major',labelsize=13)
            axPlot_g1.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r1.tick_params(axis='both',which='major',labelsize=13)
            axPlot_r1.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i1.tick_params(axis='both',which='major',labelsize=13)
            axPlot_i1.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z1.tick_params(axis='both',which='major',labelsize=13)
            axPlot_z1.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_g2.tick_params(axis='both',which='major',labelsize=13)
            axPlot_g2.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r2.tick_params(axis='both',which='major',labelsize=13)
            axPlot_r2.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i2.tick_params(axis='both',which='major',labelsize=13)
            axPlot_i2.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z2.tick_params(axis='both',which='major',labelsize=13)
            axPlot_z2.tick_params(axis='both',which='minor',labelsize=0)

            axPlot_g1.set_xticks(major_xticks)
            axPlot_g1.set_xticks(minor_xticks,minor=True)
            axPlot_r1.set_xticks(major_xticks)
            axPlot_r1.set_xticks(minor_xticks,minor=True)
            axPlot_i1.set_xticks(major_xticks)
            axPlot_i1.set_xticks(minor_xticks,minor=True)
            axPlot_z1.set_xticks(major_xticks)
            axPlot_z1.set_xticks(minor_xticks,minor=True)
            axPlot_g2.set_xticks(major_xticks)
            axPlot_g2.set_xticks(minor_xticks,minor=True)
            axPlot_r2.set_xticks(major_xticks)
            axPlot_r2.set_xticks(minor_xticks,minor=True)
            axPlot_i2.set_xticks(major_xticks)
            axPlot_i2.set_xticks(minor_xticks,minor=True)
            axPlot_z2.set_xticks(major_xticks)
            axPlot_z2.set_xticks(minor_xticks,minor=True)

            axPlot_g1.set_yticks(major_yticks)
            axPlot_g1.set_yticks(minor_yticks,minor=True)
            axPlot_r1.set_yticks(major_yticks)
            axPlot_r1.set_yticks(minor_yticks,minor=True)
            axPlot_i1.set_yticks(major_yticks)
            axPlot_i1.set_yticks(minor_yticks,minor=True)
            axPlot_z1.set_yticks(major_yticks)
            axPlot_z1.set_yticks(minor_yticks,minor=True)
            axPlot_g2.set_yticks(major_yticks)
            axPlot_g2.set_yticks(minor_yticks,minor=True)
            axPlot_r2.set_yticks(major_yticks)
            axPlot_r2.set_yticks(minor_yticks,minor=True)
            axPlot_i2.set_yticks(major_yticks)
            axPlot_i2.set_yticks(minor_yticks,minor=True)
            axPlot_z2.set_yticks(major_yticks)
            axPlot_z2.set_yticks(minor_yticks,minor=True)

            axPlot_g1.tick_params(which='both',direction='in')
            axPlot_r1.tick_params(which='both',direction='in')
            axPlot_i1.tick_params(which='both',direction='in')
            axPlot_z1.tick_params(which='both',direction='in')
            axPlot_g2.tick_params(which='both',direction='in')
            axPlot_r2.tick_params(which='both',direction='in')
            axPlot_i2.tick_params(which='both',direction='in')
            axPlot_z2.tick_params(which='both',direction='in')

            axPlot_g1.set_xlim([-0.1,0.9])
            axPlot_r1.set_xlim([-0.1,0.9])
            axPlot_i1.set_xlim([-0.1,0.9])
            axPlot_z1.set_xlim([-0.1,0.9])
            axPlot_g2.set_xlim([-0.1,0.9])
            axPlot_r2.set_xlim([-0.1,0.9])
            axPlot_i2.set_xlim([-0.1,0.9])
            axPlot_z2.set_xlim([-0.1,0.9])
            #axPlot_g.set_xlim([-1,1])
            #axPlot_r.set_xlim([-1,1])
            #axPlot_i.set_xlim([-1,1])
            #axPlot_z.set_xlim([-1,1])
            axPlot_g1.set_ylim([0,0.45])
            axPlot_r1.set_ylim([0,0.45])
            axPlot_i1.set_ylim([0,0.45])
            axPlot_z1.set_ylim([0,0.45])
            axPlot_g2.set_ylim([0,0.45])
            axPlot_r2.set_ylim([0,0.45])
            axPlot_i2.set_ylim([0,0.45])
            axPlot_z2.set_ylim([0,0.45])

            axPlot_g1.xaxis.set_ticks_position('both')
            axPlot_r1.xaxis.set_ticks_position('both')
            axPlot_i1.xaxis.set_ticks_position('both')
            axPlot_z1.xaxis.set_ticks_position('both')
            axPlot_g1.yaxis.set_ticks_position('both')
            axPlot_r1.yaxis.set_ticks_position('both')
            axPlot_i1.yaxis.set_ticks_position('both')
            axPlot_z1.yaxis.set_ticks_position('both')
            axPlot_g2.xaxis.set_ticks_position('both')
            axPlot_r2.xaxis.set_ticks_position('both')
            axPlot_i2.xaxis.set_ticks_position('both')
            axPlot_z2.xaxis.set_ticks_position('both')
            axPlot_g2.yaxis.set_ticks_position('both')
            axPlot_r2.yaxis.set_ticks_position('both')
            axPlot_i2.yaxis.set_ticks_position('both')
            axPlot_z2.yaxis.set_ticks_position('both')

            for item in p_g1:
                item.set_height(item.get_height()/sum(x_g1))
            for item in p_r1:
                item.set_height(item.get_height()/sum(x_r1))
            for item in p_i1:
                item.set_height(item.get_height()/sum(x_i1))
            for item in p_z1:
                item.set_height(item.get_height()/sum(x_z1))
            for item in p_g2:
                item.set_height(item.get_height()/sum(x_g2))
            for item in p_r2:
                item.set_height(item.get_height()/sum(x_r2))
            for item in p_i2:
                item.set_height(item.get_height()/sum(x_i2))
            for item in p_z2:
                item.set_height(item.get_height()/sum(x_z2))

            max_hist_g1,max_hist_r1,max_hist_i1,max_hist_z1 = 0,0,0,0
            for item in p_g1:
                if item.get_height() > max_hist_g1: max_hist_g1 = item.get_height()
            for item in p_r1:
                if item.get_height() > max_hist_r1: max_hist_r1 = item.get_height()
            for item in p_i1:
                if item.get_height() > max_hist_i1: max_hist_i1 = item.get_height()
            for item in p_z1:
                if item.get_height() > max_hist_z1: max_hist_z1 = item.get_height()
                    
            max_hist_g2,max_hist_r2,max_hist_i2,max_hist_z2 = 0,0,0,0
            for item in p_g2:
                if item.get_height() > max_hist_g2: max_hist_g2 = item.get_height()
            for item in p_r2:
                if item.get_height() > max_hist_r2: max_hist_r2 = item.get_height()
            for item in p_i2:
                if item.get_height() > max_hist_i2: max_hist_i2 = item.get_height()
            for item in p_z2:
                if item.get_height() > max_hist_z2: max_hist_z2 = item.get_height()
            #print(max_hist_g,max_hist_r,max_hist_i,max_hist_z)
            
            num_histbars_g1 = len(p_g1)
            num_histbars_r1 = len(p_r1)
            num_histbars_i1 = len(p_i1)
            num_histbars_z1 = len(p_z1)
            num_histbars_g2 = len(p_g2)
            num_histbars_r2 = len(p_r2)
            num_histbars_i2 = len(p_i2)
            num_histbars_z2 = len(p_z2)
            height_histbars_g1 = [0 for idx in range(num_histbars_g1)]
            height_histbars_r1 = [0 for idx in range(num_histbars_r1)]
            height_histbars_i1 = [0 for idx in range(num_histbars_i1)]
            height_histbars_z1 = [0 for idx in range(num_histbars_z1)]
            height_histbars_g2 = [0 for idx in range(num_histbars_g2)]
            height_histbars_r2 = [0 for idx in range(num_histbars_r2)]
            height_histbars_i2 = [0 for idx in range(num_histbars_i2)]
            height_histbars_z2 = [0 for idx in range(num_histbars_z2)]

            idx = 0
            for item in p_g1:
                height_histbars_g1[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_r1:
                height_histbars_r1[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_i1:
                height_histbars_i1[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_z1:
                height_histbars_z1[idx] = item.get_height()
                idx += 1
            
            idx = 0
            for item in p_g2:
                height_histbars_g2[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_r2:
                height_histbars_r2[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_i2:
                height_histbars_i2[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_z2:
                height_histbars_z2[idx] = item.get_height()
                idx += 1
            
            x_hist_g1 = [0 for idx in range(len(bins_g1)-1)]
            x_hist_r1 = [0 for idx in range(len(bins_r1)-1)]
            x_hist_i1 = [0 for idx in range(len(bins_i1)-1)]
            x_hist_z1 = [0 for idx in range(len(bins_z1)-1)]
            x_hist_g2 = [0 for idx in range(len(bins_g2)-1)]
            x_hist_r2 = [0 for idx in range(len(bins_r2)-1)]
            x_hist_i2 = [0 for idx in range(len(bins_i2)-1)]
            x_hist_z2 = [0 for idx in range(len(bins_z2)-1)]
            for idx in range(0,len(x_hist_g1)):
                x_hist_g1[idx] = (bins_g1[idx] + bins_g1[idx+1])/2
            for idx in range(0,len(x_hist_r1)):
                x_hist_r1[idx] = (bins_r1[idx] + bins_r1[idx+1])/2
            for idx in range(0,len(x_hist_i)):
                x_hist_i1[idx] = (bins_i1[idx] + bins_i1[idx+1])/2
            for idx in range(0,len(x_hist_z)):
                x_hist_z1[idx] = (bins_z1[idx] + bins_z1[idx+1])/2
            for idx in range(0,len(x_hist_g2)):
                x_hist_g2[idx] = (bins_g2[idx] + bins_g2[idx+1])/2
            for idx in range(0,len(x_hist_r2)):
                x_hist_r2[idx] = (bins_r2[idx] + bins_r2[idx+1])/2
            for idx in range(0,len(x_hist_i2)):
                x_hist_i2[idx] = (bins_i2[idx] + bins_i2[idx+1])/2
            for idx in range(0,len(x_hist_z2)):
                x_hist_z2[idx] = (bins_z2[idx] + bins_z2[idx+1])/2
            
            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g1.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g1.transAxes,bbox=props1)
            axPlot_r1.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r1.transAxes,bbox=props1)
            axPlot_i1.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i1.transAxes,bbox=props1)
            axPlot_z1.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z1.transAxes,bbox=props1)
            #axPlot_g.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_g)+'\n'+r'$\sigma={:.4f}$'.format(sigma_g),fontsize=11,transform=axPlot_g.transAxes,bbox=props2)
            #axPlot_r.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_r)+'\n'+r'$\sigma={:.4f}$'.format(sigma_r),fontsize=11,transform=axPlot_r.transAxes,bbox=props2)
            #axPlot_i.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_i)+'\n'+r'$\sigma={:.4f}$'.format(sigma_i),fontsize=11,transform=axPlot_i.transAxes,bbox=props2)
            #axPlot_z.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_z)+'\n'+r'$\sigma={:.4f}$'.format(sigma_z),fontsize=11,transform=axPlot_z.transAxes,bbox=props2)
            axPlot_g1.set_xlabel(r'$\delta m_{\rm SDSS}$',fontsize=14)
            axPlot_r1.set_xlabel(r'$\delta m_{\rm SDSS}$',fontsize=14)
            axPlot_i1.set_xlabel(r'$\delta m_{\rm SDSS}$',fontsize=14)
            axPlot_z1.set_xlabel(r'$\delta m_{\rm SDSS}$',fontsize=14)
            axPlot_g2.set_xlabel(r'$\delta m_{\rm Mcdm}$',fontsize=14)
            axPlot_r2.set_xlabel(r'$\delta m_{\rm Mcdm}$',fontsize=14)
            axPlot_i2.set_xlabel(r'$\delta m_{\rm Mcdm}$',fontsize=14)
            axPlot_z2.set_xlabel(r'$\delta m_{\rm Mcdm}$',fontsize=14)
            axPlot_g1.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_r1.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_i1.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_z1.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_g2.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_r2.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_i2.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_z2.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            
            plt.savefig(plot_filepath,format='pdf',transparent=True)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Histograms of uncertainties for MACADAMIA and SDSS data plotted successfully.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_err_histograms_griz()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_err_histograms_griz()')
    return keep_going


def create_diff_histograms_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.43,0.35,0.25]
            rect_plot_2 = [0.12,0.12,0.35,0.25]
            rect_plot_3 = [0.55,0.43,0.35,0.25]
            rect_plot_4 = [0.55,0.12,0.35,0.25]
            plt.figure(1,figsize=(12,12))
            axPlot_g = plt.axes(rect_plot_1)
            axPlot_r = plt.axes(rect_plot_2)
            axPlot_i = plt.axes(rect_plot_3)
            axPlot_z = plt.axes(rect_plot_4)
            
            n_mag_diff_bins = 81
            #mag_diff_min = np.amin(mag_diffs)
            #mag_diff_max = np.amax(mag_diffs)
            #n_mag_diff_bins = 401
            mag_diff_min = -1.0125
            mag_diff_max = 1.0125
            mag_diff_bins = np.arange(mag_diff_min,mag_diff_max,(mag_diff_max-mag_diff_min)/n_mag_diff_bins)
            mu_g,sigma_g = norm.fit(mag_diffs_g)
            mu_r,sigma_r = norm.fit(mag_diffs_r)
            mu_i,sigma_i = norm.fit(mag_diffs_i)
            mu_z,sigma_z = norm.fit(mag_diffs_z)
            x_g,bins_g,p_g = axPlot_g.hist(mag_diffs_g,bins=mag_diff_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_r,bins_r,p_r = axPlot_r.hist(mag_diffs_r,bins=mag_diff_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_i,bins_i,p_i = axPlot_i.hist(mag_diffs_i,bins=mag_diff_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            x_z,bins_z,p_z = axPlot_z.hist(mag_diffs_z,bins=mag_diff_bins,histtype='bar',density=True,edgecolor='black',color='#2077b4',zorder=2)
            
            major_xticks = np.arange(-2,2,0.2)
            minor_xticks = np.arange(-2,2,0.1)
            #major_xticks = np.arange(-12,12,2)
            #minor_xticks = np.arange(-12,12,0.5)
            major_yticks = np.arange(-2,2,0.1)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g.tick_params(axis='both',which='major',labelsize=13)
            axPlot_g.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r.tick_params(axis='both',which='major',labelsize=13)
            axPlot_r.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i.tick_params(axis='both',which='major',labelsize=13)
            axPlot_i.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z.tick_params(axis='both',which='major',labelsize=13)
            axPlot_z.tick_params(axis='both',which='minor',labelsize=0)

            axPlot_g.set_xticks(major_xticks)
            axPlot_g.set_xticks(minor_xticks,minor=True)
            axPlot_r.set_xticks(major_xticks)
            axPlot_r.set_xticks(minor_xticks,minor=True)
            axPlot_i.set_xticks(major_xticks)
            axPlot_i.set_xticks(minor_xticks,minor=True)
            axPlot_z.set_xticks(major_xticks)
            axPlot_z.set_xticks(minor_xticks,minor=True)

            axPlot_g.set_yticks(major_yticks)
            axPlot_g.set_yticks(minor_yticks,minor=True)
            axPlot_r.set_yticks(major_yticks)
            axPlot_r.set_yticks(minor_yticks,minor=True)
            axPlot_i.set_yticks(major_yticks)
            axPlot_i.set_yticks(minor_yticks,minor=True)
            axPlot_z.set_yticks(major_yticks)
            axPlot_z.set_yticks(minor_yticks,minor=True)

            axPlot_g.tick_params(which='both',direction='in')
            axPlot_r.tick_params(which='both',direction='in')
            axPlot_i.tick_params(which='both',direction='in')
            axPlot_z.tick_params(which='both',direction='in')

            axPlot_g.set_xlim([-0.5,0.5])
            axPlot_r.set_xlim([-0.5,0.5])
            axPlot_i.set_xlim([-0.5,0.5])
            axPlot_z.set_xlim([-0.5,0.5])
            #axPlot_g.set_xlim([-1,1])
            #axPlot_r.set_xlim([-1,1])
            #axPlot_i.set_xlim([-1,1])
            #axPlot_z.set_xlim([-1,1])
            axPlot_g.set_ylim([0,0.45])
            axPlot_r.set_ylim([0,0.45])
            axPlot_i.set_ylim([0,0.45])
            axPlot_z.set_ylim([0,0.45])

            axPlot_g.xaxis.set_ticks_position('both')
            axPlot_r.xaxis.set_ticks_position('both')
            axPlot_i.xaxis.set_ticks_position('both')
            axPlot_z.xaxis.set_ticks_position('both')
            axPlot_g.yaxis.set_ticks_position('both')
            axPlot_r.yaxis.set_ticks_position('both')
            axPlot_i.yaxis.set_ticks_position('both')
            axPlot_z.yaxis.set_ticks_position('both')
    
            for item in p_g:
                item.set_height(item.get_height()/sum(x_g))
            for item in p_r:
                item.set_height(item.get_height()/sum(x_r))
            for item in p_i:
                item.set_height(item.get_height()/sum(x_i))
            for item in p_z:
                item.set_height(item.get_height()/sum(x_z))

            max_hist_g,max_hist_r,max_hist_i,max_hist_z = 0,0,0,0
            for item in p_g:
                if item.get_height() > max_hist_g: max_hist_g = item.get_height()
            for item in p_r:
                if item.get_height() > max_hist_r: max_hist_r = item.get_height()
            for item in p_i:
                if item.get_height() > max_hist_i: max_hist_i = item.get_height()
            for item in p_z:
                if item.get_height() > max_hist_z: max_hist_z = item.get_height()
            #print(max_hist_g,max_hist_r,max_hist_i,max_hist_z)
            
            num_histbars_g = len(p_g)
            num_histbars_r = len(p_r)
            num_histbars_i = len(p_i)
            num_histbars_z = len(p_z)
            height_histbars_g = [0 for idx in range(num_histbars_g)]
            height_histbars_r = [0 for idx in range(num_histbars_r)]
            height_histbars_i = [0 for idx in range(num_histbars_i)]
            height_histbars_z = [0 for idx in range(num_histbars_z)]

            idx = 0
            for item in p_g:
                height_histbars_g[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_r:
                height_histbars_r[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_i:
                height_histbars_i[idx] = item.get_height()
                idx += 1
            idx = 0
            for item in p_z:
                height_histbars_z[idx] = item.get_height()
                idx += 1
            
            x_hist_g = [0 for idx in range(len(bins_g)-1)]
            x_hist_r = [0 for idx in range(len(bins_r)-1)]
            x_hist_i = [0 for idx in range(len(bins_i)-1)]
            x_hist_z = [0 for idx in range(len(bins_z)-1)]
            for idx in range(0,len(x_hist_g)):
                x_hist_g[idx] = (bins_g[idx] + bins_g[idx+1])/2
            for idx in range(0,len(x_hist_r)):
                x_hist_r[idx] = (bins_r[idx] + bins_r[idx+1])/2
            for idx in range(0,len(x_hist_i)):
                x_hist_i[idx] = (bins_i[idx] + bins_i[idx+1])/2
            for idx in range(0,len(x_hist_z)):
                x_hist_z[idx] = (bins_z[idx] + bins_z[idx+1])/2
            
            #print(x_g)
            #print(bins_g)
            #print('g')
            #for idx in range(0,len(x_hist_g)):
            #    print('{:.2f},{:.3f} '.format(x_hist_g[idx],height_histbars_g[idx]))
            #print('\nr')
            #for idx in range(0,len(x_hist_r)):
            #    print('{:.2f},{:.3f} '.format(x_hist_r[idx],height_histbars_r[idx]))
            #print('\ni')
            #for idx in range(0,len(x_hist_i)):
            #    print('{:.2f},{:.3f} '.format(x_hist_i[idx],height_histbars_i[idx]))
            #print('\nz')
            #for idx in range(0,len(x_hist_z)):
            #    print('{:.2f},{:.3f} '.format(x_hist_z[idx],height_histbars_z[idx]))
            
            popt_g,_ = optimize.curve_fit(gaussian,x_hist_g,height_histbars_g,bounds=([0.1,-0.1,0],[1,0.1,1]))
            popt_r,_ = optimize.curve_fit(gaussian,x_hist_r,height_histbars_r,bounds=([0.1,-0.1,0],[1,0.1,1]))
            popt_i,_ = optimize.curve_fit(gaussian,x_hist_i,height_histbars_i,bounds=([0.1,-0.1,0],[1,0.1,1]))
            popt_z,_ = optimize.curve_fit(gaussian,x_hist_z,height_histbars_z,bounds=([0.1,-0.1,0],[1,0.1,1]))
            
            mu_g,sigma_g = popt_g[1],popt_g[2]
            mu_r,sigma_r = popt_r[1],popt_r[2]
            mu_i,sigma_i = popt_i[1],popt_i[2]
            mu_z,sigma_z = popt_z[1],popt_z[2]
            
            x_normpdf = np.arange(-1.025,1.025,0.01)
            y_g = scipy.stats.norm.pdf(x_normpdf,mu_g,sigma_g)
            y_r = scipy.stats.norm.pdf(x_normpdf,mu_r,sigma_r)
            y_i = scipy.stats.norm.pdf(x_normpdf,mu_i,sigma_i)
            y_z = scipy.stats.norm.pdf(x_normpdf,mu_z,sigma_z)
            max_gaussian_g = max(y_g)
            max_gaussian_r = max(y_r)
            max_gaussian_i = max(y_i)
            max_gaussian_z = max(y_z)
            #print(max_gaussian_g,max_gaussian_r,max_gaussian_i,max_gaussian_z)
            
            #axPlot_g.plot(x_normpdf,y_g/(max_gaussian_g/max_hist_g),ls='--',lw=2,color='#000088',zorder=2)
            #axPlot_r.plot(x_normpdf,y_r/(max_gaussian_r/max_hist_r),ls='--',lw=2,color='#000088',zorder=2)
            #axPlot_i.plot(x_normpdf,y_i/(max_gaussian_i/max_hist_i),ls='--',lw=2,color='#000088',zorder=2)
            #axPlot_z.plot(x_normpdf,y_z/(max_gaussian_z/max_hist_z),ls='--',lw=2,color='#000088',zorder=2)
            
            axPlot_g.plot(x_normpdf,gaussian(x_normpdf,*popt_g),ls='--',lw=2,color='#880000',zorder=3)
            axPlot_r.plot(x_normpdf,gaussian(x_normpdf,*popt_r),ls='--',lw=2,color='#880000',zorder=3)
            axPlot_i.plot(x_normpdf,gaussian(x_normpdf,*popt_i),ls='--',lw=2,color='#880000',zorder=3)
            axPlot_z.plot(x_normpdf,gaussian(x_normpdf,*popt_z),ls='--',lw=2,color='#880000',zorder=3)
            #print('\nFit parameters')
            #print('{:.3f},{:.3f},{:.3f}'.format(popt_g[0],popt_g[1],popt_g[2]))
            #print('{:.3f},{:.3f},{:.3f}'.format(popt_r[0],popt_r[1],popt_r[2]))
            #print('{:.3f},{:.3f},{:.3f}'.format(popt_i[0],popt_i[1],popt_i[2]))
            #print('{:.3f},{:.3f},{:.3f}'.format(popt_z[0],popt_z[1],popt_z[2]))
            
            
            axPlot_g.axvline(x=mu_g,color='#000000',lw=1.5,ls='--',zorder=1)
            axPlot_g.axvline(x=mu_g-sigma_g,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_g.axvline(x=mu_g+sigma_g,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_g.axvline(x=mu_g-3*sigma_g,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_g.axvline(x=mu_g+3*sigma_g,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_r.axvline(x=mu_r,color='#000000',lw=1.5,ls='--',zorder=1)
            axPlot_r.axvline(x=mu_r-sigma_r,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_r.axvline(x=mu_r+sigma_r,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_r.axvline(x=mu_r-3*sigma_r,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_r.axvline(x=mu_r+3*sigma_r,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_i.axvline(x=mu_i,color='#000000',lw=1.5,ls='--',zorder=1)
            axPlot_i.axvline(x=mu_i-sigma_i,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_i.axvline(x=mu_i+sigma_i,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_i.axvline(x=mu_i-3*sigma_i,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_i.axvline(x=mu_i+3*sigma_i,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_z.axvline(x=mu_z,color='#000000',lw=1.5,ls='--',zorder=1)
            axPlot_z.axvline(x=mu_z-sigma_z,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_z.axvline(x=mu_z+sigma_z,color='#000000',lw=1.5,ls=':',zorder=1)
            axPlot_z.axvline(x=mu_z-3*sigma_z,color='#000000',lw=1,ls=':',zorder=1)
            axPlot_z.axvline(x=mu_z+3*sigma_z,color='#000000',lw=1,ls=':',zorder=1)

            #axPlot_g.set_xlim([-0.5,0.5])
            #axPlot_r.set_xlim([-0.5,0.5])
            #axPlot_i.set_xlim([-0.5,0.5])
            #axPlot_z.set_xlim([-0.5,0.5])
            #axPlot_g.set_ylim([0,0.45])
            #axPlot_r.set_ylim([0,0.45])
            #axPlot_i.set_ylim([0,0.45])
            #axPlot_z.set_ylim([0,0.45])
            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g.transAxes,bbox=props1)
            axPlot_r.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r.transAxes,bbox=props1)
            axPlot_i.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i.transAxes,bbox=props1)
            axPlot_z.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z.transAxes,bbox=props1)
            axPlot_g.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_g)+'\n'+r'$\sigma={:.4f}$'.format(sigma_g),fontsize=11,transform=axPlot_g.transAxes,bbox=props2)
            axPlot_r.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_r)+'\n'+r'$\sigma={:.4f}$'.format(sigma_r),fontsize=11,transform=axPlot_r.transAxes,bbox=props2)
            axPlot_i.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_i)+'\n'+r'$\sigma={:.4f}$'.format(sigma_i),fontsize=11,transform=axPlot_i.transAxes,bbox=props2)
            axPlot_z.text(0.045,0.80,r'$\mu={:.4f}$'.format(mu_z)+'\n'+r'$\sigma={:.4f}$'.format(sigma_z),fontsize=11,transform=axPlot_z.transAxes,bbox=props2)
            axPlot_g.set_xlabel(r'$m-m_{\rm SDSS}$',fontsize=14)
            axPlot_r.set_xlabel(r'$m-m_{\rm SDSS}$',fontsize=14)
            axPlot_i.set_xlabel(r'$m-m_{\rm SDSS}$',fontsize=14)
            axPlot_z.set_xlabel(r'$m-m_{\rm SDSS}$',fontsize=14)
            axPlot_g.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_r.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_i.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            axPlot_z.set_ylabel('Sample fraction',fontsize=14,labelpad=9)
            
            plt.savefig(plot_filepath,format='pdf',transparent=True)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude difference histogram generation and Gaussian fitting completed successfully.')
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_diff_histograms_griz()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_diff_histograms_griz()')
    return keep_going


def create_dev_vs_mag_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,mags_g,mags_r,mags_i,mags_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.43,0.35,0.25]
            rect_plot_2 = [0.12,0.12,0.35,0.25]
            rect_plot_3 = [0.55,0.43,0.35,0.25]
            rect_plot_4 = [0.55,0.12,0.35,0.25]
            plt.figure(1,figsize=(12,12))
            axPlot_g = plt.axes(rect_plot_1)
            axPlot_r = plt.axes(rect_plot_2)
            axPlot_i = plt.axes(rect_plot_3)
            axPlot_z = plt.axes(rect_plot_4)
    
            axPlot_g.scatter(mags_g,mag_diffs_g,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_r.scatter(mags_r,mag_diffs_r,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_i.scatter(mags_i,mag_diffs_i,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_z.scatter(mags_z,mag_diffs_z,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            
            major_xticks = np.arange(0,30,0.5)
            minor_xticks = np.arange(0,30,0.1)
            major_yticks = np.arange(-2,2,0.2)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g.tick_params(axis='both',which='major',labelsize=12)
            axPlot_r.tick_params(axis='both',which='major',labelsize=12)
            axPlot_i.tick_params(axis='both',which='major',labelsize=12)
            axPlot_z.tick_params(axis='both',which='major',labelsize=12)
            axPlot_g.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_g.set_xticks(major_xticks)
            axPlot_r.set_xticks(major_xticks)
            axPlot_i.set_xticks(major_xticks)
            axPlot_z.set_xticks(major_xticks)
            axPlot_g.set_xticks(minor_xticks,minor=True)
            axPlot_r.set_xticks(minor_xticks,minor=True)
            axPlot_i.set_xticks(minor_xticks,minor=True)
            axPlot_z.set_xticks(minor_xticks,minor=True)

            axPlot_g.set_yticks(major_yticks)
            axPlot_r.set_yticks(major_yticks)
            axPlot_i.set_yticks(major_yticks)
            axPlot_z.set_yticks(major_yticks)
            axPlot_g.set_yticks(minor_yticks,minor=True)
            axPlot_r.set_yticks(minor_yticks,minor=True)
            axPlot_i.set_yticks(minor_yticks,minor=True)
            axPlot_z.set_yticks(minor_yticks,minor=True)
            axPlot_g.tick_params(which='both',direction='in')
            axPlot_r.tick_params(which='both',direction='in')
            axPlot_i.tick_params(which='both',direction='in')
            axPlot_z.tick_params(which='both',direction='in')

            axPlot_g.set_xlim([16.3,20.8])
            axPlot_r.set_xlim([16.3,20.8])
            axPlot_i.set_xlim([16.3,20.8])
            axPlot_z.set_xlim([16.3,20.8])
            axPlot_g.set_ylim([-1.05,1.05])
            axPlot_r.set_ylim([-1.05,1.05])
            axPlot_i.set_ylim([-1.05,1.05])
            axPlot_z.set_ylim([-1.05,1.05])

            axPlot_g.xaxis.set_ticks_position('both')
            axPlot_r.xaxis.set_ticks_position('both')
            axPlot_i.xaxis.set_ticks_position('both')
            axPlot_z.xaxis.set_ticks_position('both')
            axPlot_g.yaxis.set_ticks_position('both')
            axPlot_r.yaxis.set_ticks_position('both')
            axPlot_i.yaxis.set_ticks_position('both')
            axPlot_z.yaxis.set_ticks_position('both')
            
            axPlot_g.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_g.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)

            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g.transAxes,bbox=props1)
            axPlot_r.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r.transAxes,bbox=props1)
            axPlot_i.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i.transAxes,bbox=props1)
            axPlot_z.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z.transAxes,bbox=props1)

            plt.savefig(plot_filepath,format='pdf',transparent=True,dpi=400)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude error vs. magnitude plot generation complete.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_dev_vs_mag_griz()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_dev_vs_mag_griz()')
    
    return keep_going

def create_dev_vs_mag_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,mags_g,mags_r,mags_i,mags_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.43,0.35,0.25]
            rect_plot_2 = [0.12,0.12,0.35,0.25]
            rect_plot_3 = [0.55,0.43,0.35,0.25]
            rect_plot_4 = [0.55,0.12,0.35,0.25]
            plt.figure(1,figsize=(12,12))
            axPlot_g = plt.axes(rect_plot_1)
            axPlot_r = plt.axes(rect_plot_2)
            axPlot_i = plt.axes(rect_plot_3)
            axPlot_z = plt.axes(rect_plot_4)
    
            axPlot_g.scatter(mags_g,mag_diffs_g,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_r.scatter(mags_r,mag_diffs_r,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_i.scatter(mags_i,mag_diffs_i,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_z.scatter(mags_z,mag_diffs_z,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            
            major_xticks = np.arange(0,30,0.5)
            minor_xticks = np.arange(0,30,0.1)
            major_yticks = np.arange(-2,2,0.2)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g.tick_params(axis='both',which='major',labelsize=12)
            axPlot_r.tick_params(axis='both',which='major',labelsize=12)
            axPlot_i.tick_params(axis='both',which='major',labelsize=12)
            axPlot_z.tick_params(axis='both',which='major',labelsize=12)
            axPlot_g.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_g.set_xticks(major_xticks)
            axPlot_r.set_xticks(major_xticks)
            axPlot_i.set_xticks(major_xticks)
            axPlot_z.set_xticks(major_xticks)
            axPlot_g.set_xticks(minor_xticks,minor=True)
            axPlot_r.set_xticks(minor_xticks,minor=True)
            axPlot_i.set_xticks(minor_xticks,minor=True)
            axPlot_z.set_xticks(minor_xticks,minor=True)

            axPlot_g.set_yticks(major_yticks)
            axPlot_r.set_yticks(major_yticks)
            axPlot_i.set_yticks(major_yticks)
            axPlot_z.set_yticks(major_yticks)
            axPlot_g.set_yticks(minor_yticks,minor=True)
            axPlot_r.set_yticks(minor_yticks,minor=True)
            axPlot_i.set_yticks(minor_yticks,minor=True)
            axPlot_z.set_yticks(minor_yticks,minor=True)
            axPlot_g.tick_params(which='both',direction='in')
            axPlot_r.tick_params(which='both',direction='in')
            axPlot_i.tick_params(which='both',direction='in')
            axPlot_z.tick_params(which='both',direction='in')

            axPlot_g.set_xlim([16.3,20.8])
            axPlot_r.set_xlim([16.3,20.8])
            axPlot_i.set_xlim([16.3,20.8])
            axPlot_z.set_xlim([16.3,20.8])
            axPlot_g.set_ylim([-1.05,1.05])
            axPlot_r.set_ylim([-1.05,1.05])
            axPlot_i.set_ylim([-1.05,1.05])
            axPlot_z.set_ylim([-1.05,1.05])

            axPlot_g.xaxis.set_ticks_position('both')
            axPlot_r.xaxis.set_ticks_position('both')
            axPlot_i.xaxis.set_ticks_position('both')
            axPlot_z.xaxis.set_ticks_position('both')
            axPlot_g.yaxis.set_ticks_position('both')
            axPlot_r.yaxis.set_ticks_position('both')
            axPlot_i.yaxis.set_ticks_position('both')
            axPlot_z.yaxis.set_ticks_position('both')
            
            axPlot_g.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_g.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)

            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g.transAxes,bbox=props1)
            axPlot_r.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r.transAxes,bbox=props1)
            axPlot_i.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i.transAxes,bbox=props1)
            axPlot_z.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z.transAxes,bbox=props1)

            plt.savefig(plot_filepath,format='pdf',transparent=True,dpi=400)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude error vs. magnitude plot generation complete.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_dev_vs_mag_griz()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_dev_vs_mag_griz()')
    
    return keep_going


def create_dev_vs_mag_griz_density(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,mags_g,mags_r,mags_i,mags_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.43,0.35,0.25]
            rect_plot_2 = [0.12,0.12,0.35,0.25]
            rect_plot_3 = [0.55,0.43,0.35,0.25]
            rect_plot_4 = [0.55,0.12,0.35,0.25]
            plt.figure(1,figsize=(12,12))
            axPlot_g = plt.axes(rect_plot_1)
            axPlot_r = plt.axes(rect_plot_2)
            axPlot_i = plt.axes(rect_plot_3)
            axPlot_z = plt.axes(rect_plot_4)
    
            axPlot_g.scatter(mags_g,mag_diffs_g,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_r.scatter(mags_r,mag_diffs_r,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_i.scatter(mags_i,mag_diffs_i,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_z.scatter(mags_z,mag_diffs_z,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            
            major_xticks = np.arange(0,30,0.5)
            minor_xticks = np.arange(0,30,0.1)
            major_yticks = np.arange(-2,2,0.2)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g.tick_params(axis='both',which='major',labelsize=12)
            axPlot_r.tick_params(axis='both',which='major',labelsize=12)
            axPlot_i.tick_params(axis='both',which='major',labelsize=12)
            axPlot_z.tick_params(axis='both',which='major',labelsize=12)
            axPlot_g.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_g.set_xticks(major_xticks)
            axPlot_r.set_xticks(major_xticks)
            axPlot_i.set_xticks(major_xticks)
            axPlot_z.set_xticks(major_xticks)
            axPlot_g.set_xticks(minor_xticks,minor=True)
            axPlot_r.set_xticks(minor_xticks,minor=True)
            axPlot_i.set_xticks(minor_xticks,minor=True)
            axPlot_z.set_xticks(minor_xticks,minor=True)

            axPlot_g.set_yticks(major_yticks)
            axPlot_r.set_yticks(major_yticks)
            axPlot_i.set_yticks(major_yticks)
            axPlot_z.set_yticks(major_yticks)
            axPlot_g.set_yticks(minor_yticks,minor=True)
            axPlot_r.set_yticks(minor_yticks,minor=True)
            axPlot_i.set_yticks(minor_yticks,minor=True)
            axPlot_z.set_yticks(minor_yticks,minor=True)
            axPlot_g.tick_params(which='both',direction='in')
            axPlot_r.tick_params(which='both',direction='in')
            axPlot_i.tick_params(which='both',direction='in')
            axPlot_z.tick_params(which='both',direction='in')

            axPlot_g.set_xlim([16.3,20.8])
            axPlot_r.set_xlim([16.3,20.8])
            axPlot_i.set_xlim([16.3,20.8])
            axPlot_z.set_xlim([16.3,20.8])
            axPlot_g.set_ylim([-1.05,1.05])
            axPlot_r.set_ylim([-1.05,1.05])
            axPlot_i.set_ylim([-1.05,1.05])
            axPlot_z.set_ylim([-1.05,1.05])

            axPlot_g.xaxis.set_ticks_position('both')
            axPlot_r.xaxis.set_ticks_position('both')
            axPlot_i.xaxis.set_ticks_position('both')
            axPlot_z.xaxis.set_ticks_position('both')
            axPlot_g.yaxis.set_ticks_position('both')
            axPlot_r.yaxis.set_ticks_position('both')
            axPlot_i.yaxis.set_ticks_position('both')
            axPlot_z.yaxis.set_ticks_position('both')
            
            axPlot_g.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_xlabel(r'$m_{\rm SDSS}$',fontsize=12)
            axPlot_g.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)

            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g.transAxes,bbox=props1)
            axPlot_r.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r.transAxes,bbox=props1)
            axPlot_i.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i.transAxes,bbox=props1)
            axPlot_z.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z.transAxes,bbox=props1)

            plt.savefig(plot_filepath,format='pdf',transparent=True,dpi=400)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude error vs. magnitude plot generation complete.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_dev_vs_mag_griz_density()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_dev_vs_mag_griz_density()')
    
    return keep_going



def create_dev_vs_snr_all(mag_diffs_all,snrs_all,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.12,0.40,0.30]
            plt.figure(1,figsize=(12,12))
            axPlot = plt.axes(rect_plot_1)

            axPlot.scatter(snrs_all,mag_diffs_all,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            #axPlot.errorbar(mags_all,mag_diffs_all,yerr=mag_diff_errs_all,fmt='.',ms=7,lw=1.5,capsize=3.5,capthick=1.5,color='#ffffff',ecolor='#000000',mec='#000000',label='Inactive, 2017',zorder=10)
    
            major_xticks = np.arange(0,150,20)
            minor_xticks = np.arange(0,150,5)
            major_yticks = np.arange(-2,2,0.2)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot.tick_params(axis='both',which='major',labelsize=12)
            axPlot.tick_params(axis='both',which='minor',labelsize=0)
            axPlot.set_xticks(major_xticks)
            axPlot.set_xticks(minor_xticks,minor=True)

            axPlot.set_yticks(major_yticks)
            axPlot.set_yticks(minor_yticks,minor=True)
            axPlot.tick_params(which='both',direction='in')

            axPlot.set_xlim([0,150])
            axPlot.set_ylim([-1.05,1.05])

            axPlot.xaxis.set_ticks_position('both')
            axPlot.yaxis.set_ticks_position('both')
            
            axPlot.set_xlabel('SNR',fontsize=12)
            axPlot.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            
            plt.savefig(plot_filepath,format='pdf',transparent=True,dpi=400)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude error vs. SNR plot generation complete.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_dev_vs_snr_all()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_dev_vs_snr_all()')
        
    return keep_going


def create_dev_vs_snr_griz(mag_diffs_g,mag_diffs_r,mag_diffs_i,mag_diffs_z,snrs_g,snrs_r,snrs_i,snrs_z,plot_filepath,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            rect_plot_1 = [0.12,0.43,0.35,0.25]
            rect_plot_2 = [0.12,0.12,0.35,0.25]
            rect_plot_3 = [0.55,0.43,0.35,0.25]
            rect_plot_4 = [0.55,0.12,0.35,0.25]
            plt.figure(1,figsize=(12,12))
            axPlot_g = plt.axes(rect_plot_1)
            axPlot_r = plt.axes(rect_plot_2)
            axPlot_i = plt.axes(rect_plot_3)
            axPlot_z = plt.axes(rect_plot_4)
    
            axPlot_g.scatter(snrs_g,mag_diffs_g,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_r.scatter(snrs_r,mag_diffs_r,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_i.scatter(snrs_i,mag_diffs_i,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            axPlot_z.scatter(snrs_z,mag_diffs_z,marker='.',color='#000000',s=5.0,rasterized=True,zorder=9)
            
            major_xticks = np.arange(0,150,20)
            minor_xticks = np.arange(0,150,5)
            major_yticks = np.arange(-2,2,0.2)
            minor_yticks = np.arange(-2,2,0.05)

            axPlot_g.tick_params(axis='both',which='major',labelsize=12)
            axPlot_r.tick_params(axis='both',which='major',labelsize=12)
            axPlot_i.tick_params(axis='both',which='major',labelsize=12)
            axPlot_z.tick_params(axis='both',which='major',labelsize=12)
            axPlot_g.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_r.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_i.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_z.tick_params(axis='both',which='minor',labelsize=0)
            axPlot_g.set_xticks(major_xticks)
            axPlot_r.set_xticks(major_xticks)
            axPlot_i.set_xticks(major_xticks)
            axPlot_z.set_xticks(major_xticks)
            axPlot_g.set_xticks(minor_xticks,minor=True)
            axPlot_r.set_xticks(minor_xticks,minor=True)
            axPlot_i.set_xticks(minor_xticks,minor=True)
            axPlot_z.set_xticks(minor_xticks,minor=True)

            axPlot_g.set_yticks(major_yticks)
            axPlot_r.set_yticks(major_yticks)
            axPlot_i.set_yticks(major_yticks)
            axPlot_z.set_yticks(major_yticks)
            axPlot_g.set_yticks(minor_yticks,minor=True)
            axPlot_r.set_yticks(minor_yticks,minor=True)
            axPlot_i.set_yticks(minor_yticks,minor=True)
            axPlot_z.set_yticks(minor_yticks,minor=True)
            axPlot_g.tick_params(which='both',direction='in')
            axPlot_r.tick_params(which='both',direction='in')
            axPlot_i.tick_params(which='both',direction='in')
            axPlot_z.tick_params(which='both',direction='in')

            axPlot_g.set_xlim([0,150])
            axPlot_r.set_xlim([0,150])
            axPlot_i.set_xlim([0,150])
            axPlot_z.set_xlim([0,150])
            axPlot_g.set_ylim([-1.05,1.05])
            axPlot_r.set_ylim([-1.05,1.05])
            axPlot_i.set_ylim([-1.05,1.05])
            axPlot_z.set_ylim([-1.05,1.05])

            axPlot_g.xaxis.set_ticks_position('both')
            axPlot_r.xaxis.set_ticks_position('both')
            axPlot_i.xaxis.set_ticks_position('both')
            axPlot_z.xaxis.set_ticks_position('both')
            axPlot_g.yaxis.set_ticks_position('both')
            axPlot_r.yaxis.set_ticks_position('both')
            axPlot_i.yaxis.set_ticks_position('both')
            axPlot_z.yaxis.set_ticks_position('both')
            
            axPlot_g.set_xlabel('SNR',fontsize=12)
            axPlot_r.set_xlabel('SNR',fontsize=12)
            axPlot_i.set_xlabel('SNR',fontsize=12)
            axPlot_z.set_xlabel('SNR',fontsize=12)
            axPlot_g.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_r.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_i.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)
            axPlot_z.set_ylabel(r'$m-m_{\rm SDSS}$',fontsize=12)

            props1 = dict(boxstyle='round',facecolor='white', alpha=1.0,edgecolor='0.00')
            props2 = dict(boxstyle='round,pad=0.05',facecolor='white',alpha=1.0,edgecolor='white')
            props3 = dict(boxstyle='round,pad=0.20',facecolor='white',alpha=1.0,edgecolor='black')
            axPlot_g.text(0.915,0.89,"g'",fontsize=16,transform=axPlot_g.transAxes,bbox=props1)
            axPlot_r.text(0.925,0.89,"r'",fontsize=16,transform=axPlot_r.transAxes,bbox=props1)
            axPlot_i.text(0.930,0.89,"i'",fontsize=16,transform=axPlot_i.transAxes,bbox=props1)
            axPlot_z.text(0.922,0.89,"z'",fontsize=16,transform=axPlot_z.transAxes,bbox=props1)

            plt.savefig(plot_filepath,format='pdf',transparent=True,dpi=400)
            plt.clf()
            plt.cla()
            plt.close()
            mcd.output_log_entry(path_logfile,'Magnitude error vs. SNR plot generation complete.')

            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: create_dev_vs_snr_griz()')
            mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: create_dev_vs_snr_griz()')
    
    return keep_going


def main():
    # Define filenames and paths
    if len(sys.argv)!=4:
        print('Usage:\n python3 pyt_VL01_compare_sdssmoc_to_macadamia.py [base_path] [sqlite_file] [output_file_date]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    output_file_date = sys.argv[3]

    keep_going = True

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('VL01_compare_sdssmoc_to_macadamia execution started','{:s} - VL01_compare_sdssmoc_to_macadamia execution started.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    
        # Create and open log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'VL01_compare_sdssmoc_to_macadamia')
        
        sdssmoc_filepath          = '/users/hhsieh/dropbox/hhsieh_db/astro_db/sdss/ADR4.dat'
        sdssmoc_output_filepath   = '/users/hhsieh/dropbox/hhsieh_db/astro_db/sdss/ADR4_numasts.txt'
        resolved_sdssmoc_name_log = '/users/hhsieh/dropbox/hhsieh_db/astro_db/sdss/ADR4_resolved_names.txt'
        if not os.path.isfile(sdssmoc_filepath):
            sdssmoc_filepath          = '/home/hhsieh/sdssmoc/ADR4.dat'
            sdssmoc_output_filepath   = '/home/hhsieh/sdssmoc/ADR4_numasts.txt'
            resolved_sdssmoc_name_log = '/home/hhsieh/sdssmoc/ADR4_resolved_names.txt'
        if not os.path.isfile(sdssmoc_filepath):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'SDSS MOC file {:s} not found.'.format(sdssmoc_filepath))
            mcd.send_status_email('VL01_compare_sdssmoc_to_macadamia execution failed','{:s} - VL01_compare_sdssmoc_to_macadamia execution failed - SDSS MOC file {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),sdssmoc_filepath))
            keep_going = False

    if keep_going:
        
        macadamia_output_filepath  = base_path + 'macadamia_archive_{:s}.dat.gz'.format(output_file_date)
        sdss_vs_macadamia_filepath = base_path + 'sdss_vs_macadamia_{:s}.dat'.format(output_file_date)       
        #macadamia_output_filepath  = base_path + 'macadamia_archive_20201027_first1000lines.dat.gz'
        #sdss_vs_macadamia_filepath = base_path + 'sdss_vs_macadamia_20201027_first1000lines.dat'
            
        #keep_going = parse_sdssmoc_table(sdssmoc_filepath,sdssmoc_output_filepath,resolved_sdssmoc_name_log,keep_going,path_logfile,path_errorfile)
        sdssmoc_table,keep_going = create_sdssmoc_table(sdssmoc_output_filepath,keep_going,path_logfile,path_errorfile)
        keep_going = create_sdss_macadamia_comparison_files2(sdssmoc_table,macadamia_output_filepath,sdss_vs_macadamia_filepath,keep_going,path_logfile,path_errorfile)            
            
        keep_going = create_sdss_macadamia_comparison_files(sdssmoc_table,macadamia_output_filepath,sdss_vs_macadamia_filepath,keep_going,path_logfile,path_errorfile)
        keep_going = create_sdss_macadamia_comparison_plots(sdss_vs_macadamia_filepath,'full_dataset',keep_going,path_logfile,path_errorfile)
        keep_going = create_sdss_macadamia_comparison_plots(sdss_vs_macadamia_filepath,'outliers_removed',keep_going,path_logfile,path_errorfile)
        keep_going = create_sdss_macadamia_comparison_plots(sdss_vs_macadamia_filepath,'high_precision',keep_going,path_logfile,path_errorfile)
            
        keep_going = create_sdss_macadamia_comparison_plots_outliers_removed(sdss_vs_macadamia_filepath,keep_going,path_logfile,path_errorfile)

    mcd.send_status_email('VL01_compare_sdssmoc_to_macadamia execution complete.','{:s} - VL01_compare_sdssmoc_to_macadamia execution complete.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()


