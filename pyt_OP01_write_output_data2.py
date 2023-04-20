import sys
import glob, os, subprocess
import os.path
import datetime
import sqlite3
from sqlite3 import Error
import macadamia_functions as mcd

# Edited 2021-11-07: changed code and error logging functionality
# Edited 2021-11-12: added code to replace spaces with underscores in telescope and instrument names in output
# Created 2021-11-13: output separate files for different filters to speed up lookup task later

##### FUNCTION DEFINITIONS -- DATA RETRIEVAL #####

def str_replace_spaces(string_orig):
    string_nospaces = ''
    for idx in range(0,len(string_orig)):
        if string_orig[idx] == ' ': string_nospaces = string_nospaces + '_'
        else:                       string_nospaces = string_nospaces + string_orig[idx]
    return string_nospaces

def retrieve_output_data(output_data_filepath,sqlite_file,path_logfile,path_errorfile):
    try:
        output_data_g_filepath                  = output_data_filepath[:-4] + '_g.dat'
        output_data_r_filepath                  = output_data_filepath[:-4] + '_r.dat'
        output_data_i_filepath                  = output_data_filepath[:-4] + '_i.dat'
        output_data_z_filepath                  = output_data_filepath[:-4] + '_z.dat'
        output_data_detections_filepath         = output_data_filepath[:-4] + '_detections.dat'
        output_data_detections_g_filepath       = output_data_filepath[:-4] + '_detections_g.dat'
        output_data_detections_r_filepath       = output_data_filepath[:-4] + '_detections_r.dat'
        output_data_detections_i_filepath       = output_data_filepath[:-4] + '_detections_i.dat'
        output_data_detections_z_filepath       = output_data_filepath[:-4] + '_detections_z.dat'
        output_data_nondetections_filepath      = output_data_filepath[:-4] + '_nondetections.dat'
        output_data_nondetections_g_filepath    = output_data_filepath[:-4] + '_nondetections_g.dat'
        output_data_nondetections_r_filepath    = output_data_filepath[:-4] + '_nondetections_r.dat'
        output_data_nondetections_i_filepath    = output_data_filepath[:-4] + '_nondetections_i.dat'
        output_data_nondetections_z_filepath    = output_data_filepath[:-4] + '_nondetections_z.dat'
        #output_data_filepath_gz               = output_data_filepath + '.gz'
        #output_data_detections_filepath_gz    = output_data_detections_filepath + '.gz'
        #output_data_nondetections_filepath_gz = output_data_nondetections_filepath + '.gz'
        #if os.path.exists(output_data_filepath_gz):               os.remove(output_data_filepath_gz)
        #if os.path.exists(output_data_detections_filepath_gz):    os.remove(output_data_detections_filepath_gz)
        #if os.path.exists(output_data_nondetections_filepath_gz): os.remove(output_data_nondetections_filepath_gz)
        header_line = 'detection_id   desig_num   telescope                        instrument                   \
tel_apsize   obs_code   image_ext   filter       obs_date     \
expstart_tai      expstart_jd   exptime   airmass    \
trackrate_ra   trackrate_dec    zpoint   zpoint_err   zpt_nstars   \
limit_mag_ps   limit_mag_sb   psf_width_mean  src_density   \
wcs_nstars   wcs_pixscale_x   wcs_pixscale_y   xcoord_pred   ycoord_pred     \
xcoord     ycoord   right_ascension     declination   dist_pred_postn       \
flux      dflux       mag   mag_err   psf_fwhm        SNR   \
trail_length   trail_fwhm   trail_phi       \
ra_src_1      dec_src_1    dist_src_pred_1   dist_src_actual_1   mag_src_1   magerr_src_1       \
ra_src_2      dec_src_2    dist_src_pred_2   dist_src_actual_2   mag_src_2   magerr_src_2       \
ra_src_3      dec_src_3    dist_src_pred_3   dist_src_actual_3   mag_src_3   magerr_src_3   \
local_srcdn   dist_edge_l   dist_edge_r   dist_edge_b   dist_edge_t   \
ra_predicted   dec_predicted   ra_predicted_err   dec_predicted_err         ra_rate        dec_rate   trail_pa_pred   trail_len_pred   \
eclip_long   eclip_lat        helio_dist          geo_dist   solar_elong    \
phs_ang   lunar_elong   lunar_illum   pa_antisolar   pa_neg_hel_v   orb_pl_angle   \
gal_long    gal_lat   true_anom   vmag_pred   \
flux_t_05rK   dflux_t_05rK   mag_t_05rK   magerr_t_05rK   \
flux_t_10rK   dflux_t_10rK   mag_t_10rK   magerr_t_10rK   \
flux_t_20rK   dflux_t_20rK   mag_t_20rK   magerr_t_20rK   \
flux_t_30rK   dflux_t_30rK   mag_t_30rK   magerr_t_30rK   \
flux_t_40rK   dflux_t_40rK   mag_t_40rK   magerr_t_40rK   \
flux_t_50rK   dflux_t_50rK   mag_t_50rK   magerr_t_50rK    \
h_mag   g_param   epoch_mjd       semimaj    eccentr     inclntn     argperi     \
ascnode      meananm   raw_data_link                                                                                  \
proc_data_file                          preview_img_file\n'
        with open(output_data_filepath,'w') as of:                 of.write(header_line)
        with open(output_data_g_filepath,'w') as of:               of.write(header_line)
        with open(output_data_r_filepath,'w') as of:               of.write(header_line)
        with open(output_data_i_filepath,'w') as of:               of.write(header_line)
        with open(output_data_z_filepath,'w') as of:               of.write(header_line)
        with open(output_data_detections_filepath,'w') as of:      of.write(header_line)
        with open(output_data_detections_g_filepath,'w') as of:    of.write(header_line)
        with open(output_data_detections_r_filepath,'w') as of:    of.write(header_line)
        with open(output_data_detections_i_filepath,'w') as of:    of.write(header_line)
        with open(output_data_detections_z_filepath,'w') as of:    of.write(header_line)
        with open(output_data_nondetections_filepath,'w') as of:   of.write(header_line)
        with open(output_data_nondetections_g_filepath,'w') as of: of.write(header_line)
        with open(output_data_nondetections_r_filepath,'w') as of: of.write(header_line)
        with open(output_data_nondetections_i_filepath,'w') as of: of.write(header_line)
        with open(output_data_nondetections_z_filepath,'w') as of: of.write(header_line)
        conn = mcd.create_connection(sqlite_file)
        cursor = conn.cursor()
        query = "SELECT det.detection_id,sso.desig_number,tel.telescope_name,inst.instrument_name,\
tel.aperture_size,tel.observatory_code,me.mosaic_element_num,exp.filter_name,exp.date_tai,\
exp.exposure_start_tai,exp.exposure_start_jd,exp.exposure_time,exp.airmass,\
exp.tracking_rate_ra,exp.tracking_rate_dec,exp.zero_point,exp.zpoint_err,exp.zpoint_nstars,\
exp.limiting_mag_ps_3sig,exp.limiting_mag_sb_3sig,exp.psf_width_mean,exp.source_density,\
exp.wcs_nstars,exp.wcs_pixscale_x,exp.wcs_pixscale_y,ddat.x_pixcoord_predicted,y_pixcoord_predicted,\
ddat.x_pixcoord,ddat.y_pixcoord,ddat.right_ascension,ddat.declination,ddat.dist_predicted_position,\
ddat.flux,ddat.dflux,ddat.mag,ddat.mag_err,ddat.psf_fwhm,ddat.signal_to_noise,\
ddat.trail_length,ddat.trail_fwhm,ddat.trail_phi,\
ddat.ra_source_1,ddat.dec_source_1,ddat.dist_source_pred_1,ddat.dist_source_actual_1,ddat.mag_source_1,ddat.magerr_source_1,\
ddat.ra_source_2,ddat.dec_source_2,ddat.dist_source_pred_2,ddat.dist_source_actual_2,ddat.mag_source_2,ddat.magerr_source_2,\
ddat.ra_source_3,ddat.dec_source_3,ddat.dist_source_pred_3,ddat.dist_source_actual_3,ddat.mag_source_3,ddat.magerr_source_3,\
ddat.source_density_1arcmin,ddat.dist_edge_left,ddat.dist_edge_right,ddat.dist_edge_bottom,ddat.dist_edge_top,\
ddet.ra_predicted,ddet.dec_predicted,ddet.ra_predicted_err,ddet.dec_predicted_err,ddet.ra_rate,ddet.dec_rate,ddet.trail_pa_expected,ddet.trail_length_expected,\
ddet.ecliptic_longitude,ddet.ecliptic_latitude,ddet.heliocentric_dist,ddet.geocentric_dist,ddet.solar_elongation,\
ddet.phase_angle,ddet.lunar_elongation,ddet.lunar_illumination,ddet.pa_antisolar,ddet.pa_neg_heliocentric_vel,\
ddet.orb_plane_angle,ddet.galactic_longitude,ddet.galactic_latitude,ddet.true_anomaly,ddet.vmag_predicted,\
mult.flux_t_05rKron,mult.dflux_t_05rKron,mult.mag_t_05rKron,mult.magerr_t_05rKron,\
mult.flux_t_10rKron,mult.dflux_t_10rKron,mult.mag_t_10rKron,mult.magerr_t_10rKron,\
mult.flux_t_20rKron,mult.dflux_t_20rKron,mult.mag_t_20rKron,mult.magerr_t_20rKron,\
mult.flux_t_30rKron,mult.dflux_t_30rKron,mult.mag_t_30rKron,mult.magerr_t_30rKron,\
mult.flux_t_40rKron,mult.dflux_t_40rKron,mult.mag_t_40rKron,mult.magerr_t_40rKron,\
mult.flux_t_50rKron,mult.dflux_t_50rKron,mult.mag_t_50rKron,mult.magerr_t_50rKron,\
sso.h_mag,sso.g_param,sso.epoch_mjd,sso.semimaj_axis,sso.eccentricity,sso.inclination,sso.arg_perihelion,\
sso.long_ascnode,sso.mean_anomaly,\
exp.raw_data_link,exp.proc_data_file,ddat.preview_image_file \
FROM detections AS det \
INNER JOIN detection_data AS ddat ON det.detection_id = ddat.detection_id \
INNER JOIN detection_details AS ddet ON det.detection_id = ddet.detection_id \
INNER JOIN detection_multiap_photometry AS mult ON det.detection_id = mult.detection_id \
INNER JOIN exposures AS exp ON ddat.exposure_id = exp.exposure_id \
INNER JOIN ssobjects AS sso ON ddat.ssobject_id = sso.ssobject_id \
INNER JOIN mosaic_elements AS me ON exp.mosaic_element_id = me.mosaic_element_id \
INNER JOIN instruments AS inst ON me.instrument_id = inst.instrument_id \
INNER JOIN telescopes AS tel ON inst.telescope_id = tel.telescope_id \
WHERE ddat.exposure_id != -1"
        mcd.output_log_entry(path_logfile,query)
        cursor.execute(query)
        rows = cursor.fetchall()
        if rows != None:
            with open(output_data_filepath,'a') as of, open(output_data_r_filepath,'a') as of_r, open(output_data_i_filepath,'a') as of_i, open(output_data_z_filepath,'a') as of_z:
                for row in rows:
                    det_id            = int(row[0])
                    desig_number      = int(row[1])
                    telescope         = str_replace_spaces(row[2])
                    instrument        = str_replace_spaces(row[3])
                    tel_apsize        = float(row[4])
                    obs_code          = row[5]
                    image_ext         = int(row[6])
                    filter_name       = row[7]
                    obs_date          = row[8]
                    expstart_tai      = row[9]
                    expstart_jd       = float(row[10])
                    exptime           = float(row[11])
                    airmass           = float(row[12])
                    trackrate_ra      = float(row[13])
                    trackrate_dec     = float(row[14])
                    zpoint            = float(row[15])
                    zpoint_err        = float(row[16])
                    zpt_nstars        = int(row[17])
                    limit_mag_ps      = float(row[18])
                    limit_mag_sb      = float(row[19])
                    psf_width_mean    = float(row[20])
                    src_density       = float(row[21])
                    wcs_nstars        = int(row[22])
                    wcs_pixscale_x    = float(row[23])
                    wcs_pixscale_y    = float(row[24])
                    xcoord_pred       = float(row[25])
                    ycoord_pred       = float(row[26])
                    xcoord            = float(row[27])
                    ycoord            = float(row[28])
                    right_ascension   = float(row[29])
                    declination       = float(row[30])
                    dist_pred_postn   = float(row[31])
                    flux              = float(row[32])
                    dflux             = float(row[33])
                    mag               = float(row[34])
                    mag_err           = float(row[35])
                    psf_fwhm          = float(row[36])
                    SNR               = float(row[37])
                    trail_length      = float(row[38])
                    trail_fwhm        = float(row[39])
                    trail_phi         = float(row[40])
                    ra_src_1          = float(row[41])
                    dec_src_1         = float(row[42])
                    dist_src_pred_1   = float(row[43])
                    dist_src_actual_1 = float(row[44])
                    mag_src_1         = float(row[45])
                    magerr_src_1      = float(row[46])
                    ra_src_2          = float(row[47])
                    dec_src_2         = float(row[48])
                    dist_src_pred_2   = float(row[49])
                    dist_src_actual_2 = float(row[50])
                    mag_src_2         = float(row[51])
                    magerr_src_2      = float(row[52])
                    ra_src_3          = float(row[53])
                    dec_src_3         = float(row[54])
                    dist_src_pred_3   = float(row[55])
                    dist_src_actual_3 = float(row[56])
                    mag_src_3         = float(row[57])
                    magerr_src_3      = float(row[58])
                    src_density       = float(row[59])
                    dist_edge_l       = float(row[60])
                    dist_edge_r       = float(row[61])
                    dist_edge_b       = float(row[62])
                    dist_edge_t       = float(row[63])
                    ra_predicted      = float(row[64])
                    dec_predicted     = float(row[65])
                    ra_predicted_err  = float(row[66])
                    dec_predicted_err = float(row[67])
                    ra_rate           = float(row[68])
                    dec_rate          = float(row[69])
                    trail_pa_pred     = float(row[70])
                    trail_len_pred    = float(row[71])
                    eclip_long        = float(row[72])
                    eclip_lat         = float(row[73])
                    helio_dist        = float(row[74])
                    geo_dist          = float(row[75])
                    solar_elong       = float(row[76])
                    phs_ang           = float(row[77])
                    lunar_elong       = float(row[78])
                    lunar_illum       = float(row[79])
                    pa_antisolar      = float(row[80])
                    pa_neg_hel_v      = float(row[81])
                    orbit_pl_angle    = float(row[82])
                    gal_long          = float(row[83])
                    gal_lat           = float(row[84])
                    true_anom         = float(row[85])
                    vmag_pred         = float(row[86])
                    flux_t_05rK       = float(row[87])
                    dflux_t_05rK      = float(row[88])
                    mag_t_05rK        = float(row[89])
                    magerr_t_05rK     = float(row[90])
                    flux_t_10rK       = float(row[91])
                    dflux_t_10rK      = float(row[92])
                    mag_t_10rK        = float(row[93])
                    magerr_t_10rK     = float(row[94])
                    flux_t_20rK       = float(row[95])
                    dflux_t_20rK      = float(row[96])
                    mag_t_20rK        = float(row[97])
                    magerr_t_20rK     = float(row[98])
                    flux_t_30rK       = float(row[99])
                    dflux_t_30rK      = float(row[100])
                    mag_t_30rK        = float(row[101])
                    magerr_t_30rK     = float(row[102])
                    flux_t_40rK       = float(row[103])
                    dflux_t_40rK      = float(row[104])
                    mag_t_40rK        = float(row[105])
                    magerr_t_40rK     = float(row[106])
                    flux_t_50rK       = float(row[107])
                    dflux_t_50rK      = float(row[108])
                    mag_t_50rK        = float(row[109])
                    magerr_t_50rK     = float(row[110])
                    h_mag             = float(row[111])
                    g_param           = float(row[112])
                    epoch_mjd         = float(row[113])
                    semimaj           = float(row[114])
                    eccentr           = float(row[115])
                    inclntn           = float(row[116])
                    argperi           = float(row[117])
                    ascnode           = float(row[118])
                    meananm           = float(row[119])
                    raw_data_link     = row[120]
                    proc_data_file    = row[121]
                    preview_img_file  = row[122]
                    output_line = '{:>12d}     {:>7d}   {:<30s}   {:<30s}   \
{:6.2f}        {:<4s}        {:>3d}   {:<10s}   {:10s}   \
{:>12s}   {:14.6f}    {:6.1f}     {:5.3f}   \
{:13.7f}   {:13.7f}    {:6.3f}       {:6.3f}          {:>3d}         \
{:6.3f}         {:6.3f}           {:6.3f}      {:7.3f}          \
{:>3d}          {:7.4f}          {:7.4f}      {:8.2f}      {:8.2f}   \
{:8.2f}   {:8.2f}     {:13.8f}   {:13.8f}          {:8.3f}   \
{:8.1f}   {:8.1f}  {:8.3f}  {:8.3f}   {:8.3f}   {:8.2f}       \
{:8.3f}     {:8.3f}     {:7.2f}   \
{:12.8f}   {:12.8f}           {:8.3f}            {:8.3f}    {:8.3f}       {:8.3f}   \
{:12.8f}   {:12.8f}           {:8.3f}            {:8.3f}    {:8.3f}       {:8.3f}   \
{:12.8f}   {:12.8f}           {:8.3f}            {:8.3f}    {:8.3f}       {:8.3f}       \
{:7.2f}       {:7.1f}       {:7.1f}       {:7.1f}       {:7.1f}   \
{:12.8f}    {:12.8f}       {:12.8f}        {:12.8f}   {:13.7f}   {:13.7f}         {:7.3f}          {:7.3f}     \
{:8.4f}    {:8.4f}   {:15.8f}   {:15.8f}      {:8.4f}   \
{:8.4f}        {:6.2f}         {:5.2f}        {:7.3f}        {:7.3f}        {:7.3f}   \
{:8.4f}   {:8.4f}    {:8.4f}       {:5.2f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}      \
{:8.1f}       {:8.1f}     {:8.3f}        {:8.3f}   \
{:6.3f}    {:6.3f}     {:7.1f}   {:11.6f}   {:8.6f}   {:9.5f}   {:9.5f}   \
{:9.5f}   {:10.6f}   \
{:s}   {:s}   {:s}\n'.format(det_id,desig_number,telescope,instrument,\
tel_apsize,obs_code,image_ext,filter_name,obs_date,\
expstart_tai,expstart_jd,exptime,airmass,\
trackrate_ra,trackrate_dec,zpoint,zpoint_err,zpt_nstars,\
limit_mag_ps,limit_mag_sb,psf_width_mean,src_density,\
wcs_nstars,wcs_pixscale_x,wcs_pixscale_y,xcoord_pred,ycoord_pred,\
xcoord,ycoord,right_ascension,declination,dist_pred_postn,   \
flux,dflux,mag,mag_err,psf_fwhm,SNR,   \
trail_length,trail_fwhm,trail_phi,   \
ra_src_1,dec_src_1,dist_src_pred_1,dist_src_actual_1,mag_src_1,magerr_src_1,   \
ra_src_2,dec_src_2,dist_src_pred_2,dist_src_actual_2,mag_src_2,magerr_src_2,   \
ra_src_3,dec_src_3,dist_src_pred_3,dist_src_actual_3,mag_src_3,magerr_src_3,   \
src_density,dist_edge_l,dist_edge_r,dist_edge_b,dist_edge_t,   \
ra_predicted,dec_predicted,ra_predicted_err,dec_predicted_err,ra_rate,dec_rate,trail_pa_pred,trail_len_pred,   \
eclip_long,eclip_lat,helio_dist,geo_dist,solar_elong,   \
phs_ang,lunar_elong,lunar_illum,pa_antisolar,pa_neg_hel_v,orbit_pl_angle,   \
gal_long,gal_lat,true_anom,vmag_pred,   \
flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK,\
flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK,\
flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK,\
flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK,\
flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK,\
flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK,\
h_mag,g_param,epoch_mjd,semimaj,eccentr,inclntn,argperi,\
ascnode,meananm,\
raw_data_link,proc_data_file,preview_img_file)

                    # write all output to general output files
                    of.write(output_line)
                    if mag != 99.99:
                        with open(output_data_detections_filepath,'a') as of_d: of_d.write(output_line)
                    else:
                        with open(output_data_nondetections_filepath,'a') as of_nd: of_nd.write(output_line)
                    # write g-band output to g-band output files
                    if filter_name == 'g':
                        with open(output_data_g_filepath,'a') as of_g: of_g.write(output_line)
                        if mag != 99.99:
                            with open(output_data_detections_g_filepath,'a') as of_g_d: of_g_d.write(output_line)
                        else:
                            with open(output_data_nondetections_g_filepath,'a') as of_g_nd: of_g_nd.write(output_line)
                    # write r-band output to r-band output files
                    if filter_name == 'r':
                        with open(output_data_r_filepath,'a') as of_r: of_r.write(output_line)
                        if mag != 99.99:
                            with open(output_data_detections_r_filepath,'a') as of_r_d: of_r_d.write(output_line)
                        else:
                            with open(output_data_nondetections_r_filepath,'a') as of_r_nd: of_r_nd.write(output_line)
                    # write i-band output to i-band output files
                    if filter_name == 'i':
                        with open(output_data_i_filepath,'a') as of_i: of_i.write(output_line)
                        if mag != 99.99:
                            with open(output_data_detections_i_filepath,'a') as of_i_d: of_i_d.write(output_line)
                        else:
                            with open(output_data_nondetections_i_filepath,'a') as of_i_nd: of_i_nd.write(output_line)
                    # write z-band output to z-band output files
                    if filter_name == 'z':
                        with open(output_data_z_filepath,'a') as of_z: of_z.write(output_line)
                        if mag != 99.99:
                            with open(output_data_detections_z_filepath,'a') as of_z_d: of_z_d.write(output_line)
                        else:
                            with open(output_data_nondetections_z_filepath,'a') as of_z_nd: of_z_nd.write(output_line)
        conn.close()  # Close connection to database file
    except Exception as e:
        mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
    return None


def main():
    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_ingest_detection_exposure_object_links.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('OP01_write_output_data execution started','{:s} - OP01_write_output_data execution started.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
        # Create and open log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'OP01_write_output_data')

        os.chdir(base_path)
        output_data_filepath = base_path + 'macadamia_archive_{:s}.dat'.format(datetime.datetime.today().strftime('%Y%m%d'))
        mcd.output_log_entry(path_logfile,'Writing all photometric data to {:s} ...'.format(output_data_filepath))
        retrieve_output_data(output_data_filepath,sqlite_file,path_logfile,path_errorfile)
        mcd.output_log_entry(path_logfile,'Writing photometric data to output file complete.')
        #mcd.compress_file_gzip(output_data_filepath)
        
    mcd.send_status_email('OP01_write_output_data execution complete','{:s} - OP01_write_output_data execution complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None


if __name__ == '__main__':
    main()


