import sys
import glob, os, bz2, subprocess
import os.path
import datetime
import sqlite3
from sqlite3 import Error
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE FUNCTIONS #####

def ingest_detection_data_file(detection_data_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile):
    try:
        with open(detection_data_file_to_ingest) as input_file:
            conn = mcd.create_connection(sqlite_file)
            cursor = conn.cursor()
            for _ in range(1): #skip one header line
                next(input_file)
            for line in input_file:
                detection_data  = line.split()
                detection_id            = int(detection_data[0])
                x_pixcoord_predicted    = float(detection_data[1])
                y_pixcoord_predicted    = float(detection_data[2])
                x_pixcoord              = float(detection_data[3])
                y_pixcoord              = float(detection_data[4])
                right_ascension         = float(detection_data[5])
                declination             = float(detection_data[6])
                dist_predicted_position = float(detection_data[7])
                flux_s                  = float(detection_data[8])
                dflux_s                 = float(detection_data[9])
                flux_t                  = float(detection_data[10])
                dflux_t                 = float(detection_data[11])
                mag_s                   = float(detection_data[12])
                mag_err_s               = float(detection_data[13])
                mag_t                   = float(detection_data[14])
                mag_err_t               = float(detection_data[15])
                psf_fwhm                = float(detection_data[16])
                signal_to_noise         = float(detection_data[17])
                trail_length            = float(detection_data[18])
                trail_fwhm              = float(detection_data[19])
                trail_phi               = float(detection_data[20])
                flux_t_05rK             = float(detection_data[21])
                dflux_t_05rK            = float(detection_data[22])
                mag_t_05rK              = float(detection_data[23])
                magerr_t_05rK           = float(detection_data[24])
                flux_t_10rK             = float(detection_data[25])
                dflux_t_10rK            = float(detection_data[26])
                mag_t_10rK              = float(detection_data[27])
                magerr_t_10rK           = float(detection_data[28])
                flux_t_20rK             = float(detection_data[29])
                dflux_t_20rK            = float(detection_data[30])
                mag_t_20rK              = float(detection_data[31])
                magerr_t_20rK           = float(detection_data[32])
                flux_t_30rK             = float(detection_data[33])
                dflux_t_30rK            = float(detection_data[34])
                mag_t_30rK              = float(detection_data[35])
                magerr_t_30rK           = float(detection_data[36])
                flux_t_40rK             = float(detection_data[37])
                dflux_t_40rK            = float(detection_data[38])
                mag_t_40rK              = float(detection_data[39])
                magerr_t_40rK           = float(detection_data[40])
                flux_t_50rK             = float(detection_data[41])
                dflux_t_50rK            = float(detection_data[42])
                mag_t_50rK              = float(detection_data[43])
                magerr_t_50rK           = float(detection_data[44])
                ra_source_1             = float(detection_data[45])
                dec_source_1            = float(detection_data[46])
                dist_source_pred_1      = float(detection_data[47])
                dist_source_actual_1    = float(detection_data[48])
                mag_source_1            = float(detection_data[49])
                magerr_source_1         = float(detection_data[50])
                ra_source_2             = float(detection_data[51])
                dec_source_2            = float(detection_data[52])
                dist_source_pred_2      = float(detection_data[53])
                dist_source_actual_2    = float(detection_data[54])
                mag_source_2            = float(detection_data[55])
                magerr_source_2         = float(detection_data[56])
                ra_source_3             = float(detection_data[57])
                dec_source_3            = float(detection_data[58])
                dist_source_pred_3      = float(detection_data[59])
                dist_source_actual_3    = float(detection_data[60])
                mag_source_3            = float(detection_data[61])
                magerr_source_3         = float(detection_data[62])
                source_density_1arcmin  = float(detection_data[63])
                dist_edge_left          = float(detection_data[64])
                dist_edge_right         = float(detection_data[65])
                dist_edge_bottom        = float(detection_data[66])
                dist_edge_top           = float(detection_data[67])
                preview_image_file      = detection_data[68]
                preview_image_path      = detection_data[69]
                if x_pixcoord == -999.0 and y_pixcoord == -999.0:
                    status = 'No source at predicted object position.'
                else:
                    status = 'Photometry complete.'

                mcd.output_log_entry(path_logfile,'Inserting/updating detection_data entry for detection_id {:d}...'.format(detection_id))
                # insert mags from 2D waussian fits
                query = "UPDATE detection_data SET x_pixcoord_predicted={:.3f},y_pixcoord_predicted={:.3f},x_pixcoord={:.3f},y_pixcoord={:.3f},right_ascension={:.8f},declination={:.8f},dist_predicted_position={:.3f},flux={:.6f},dflux={:.6f},mag={:.3f},mag_err={:.3f},psf_fwhm={:.3f},signal_to_noise={:.3f},trail_length={:.3f},trail_fwhm={:.3f},trail_phi={:.3f},ra_source_1={:.8f},dec_source_1={:.8f},dist_source_pred_1={:.3f},dist_source_actual_1={:.3f},mag_source_1={:.3f},magerr_source_1={:.3f},ra_source_2={:.8f},dec_source_2={:.8f},dist_source_pred_2={:.3f},dist_source_actual_2={:.3f},mag_source_2={:.3f},magerr_source_2={:.3f},ra_source_3={:.8f},dec_source_3={:.8f},dist_source_pred_3={:.3f},dist_source_actual_3={:.3f},mag_source_3={:.3f},magerr_source_3={:.3f},source_density_1arcmin={:.3f},dist_edge_left={:.1f},dist_edge_right={:.1f},dist_edge_bottom={:.1f},dist_edge_top={:.1f},preview_image_path='{:s}',preview_image_file='{:s}',detection_status='{:s}' WHERE detection_id={:d}".format(x_pixcoord_predicted,y_pixcoord_predicted,x_pixcoord,y_pixcoord,right_ascension,declination,dist_predicted_position,flux_s,dflux_s,mag_s,mag_err_s,psf_fwhm,signal_to_noise,trail_length,trail_fwhm,trail_phi,ra_source_1,dec_source_1,dist_source_pred_1,dist_source_actual_1,mag_source_1,magerr_source_1,ra_source_2,dec_source_2,dist_source_pred_2,dist_source_actual_2,mag_source_2,magerr_source_2,ra_source_3,dec_source_3,dist_source_pred_3,dist_source_actual_3,mag_source_3,magerr_source_3,source_density_1arcmin,dist_edge_left,dist_edge_right,dist_edge_bottom,dist_edge_top,preview_image_path,preview_image_file,status,detection_id)
                # insert mags from trailed waussian fits
                #query = "UPDATE detection_data SET x_pixcoord_predicted={:.3f},y_pixcoord_predicted={:.3f},x_pixcoord={:.3f},y_pixcoord={:.3f},right_ascension={:.8f},declination={:.8f},dist_predicted_position={:.3f},flux={:.6f},dflux={:.6f},mag={:.3f},mag_err={:.3f},psf_fwhm={:.3f},signal_to_noise={:.3f},trail_length={:.3f},trail_fwhm={:.3f},trail_phi={:.3f},ra_source_1={:.8f},dec_source_1={:.8f},dist_source_pred_1={:.3f},dist_source_actual_1={:.3f},mag_source_1={:.3f},magerr_source_1={:.3f},ra_source_2={:.8f},dec_source_2={:.8f},dist_source_pred_2={:.3f},dist_source_actual_2={:.3f},mag_source_2={:.3f},magerr_source_2={:.3f},ra_source_3={:.8f},dec_source_3={:.8f},dist_source_pred_3={:.3f},dist_source_actual_3={:.3f},mag_source_3={:.3f},magerr_source_3={:.3f},source_density_1arcmin={:.3f},dist_edge_left={:.1f},dist_edge_right={:.1f},dist_edge_bottom={:.1f},dist_edge_top={:.1f},preview_image_path='{:s}',preview_image_file='{:s}',detection_status='{:s}' WHERE detection_id={:d}".format(x_pixcoord_predicted,y_pixcoord_predicted,x_pixcoord,y_pixcoord,right_ascension,declination,dist_predicted_position,flux_t,dflux_t,mag_t,mag_err_t,psf_fwhm,signal_to_noise,trail_length,trail_fwhm,trail_phi,ra_source_1,dec_source_1,dist_source_pred_1,dist_source_actual_1,mag_source_1,magerr_source_1,ra_source_2,dec_source_2,dist_source_pred_2,dist_source_actual_2,mag_source_2,magerr_source_2,ra_source_3,dec_source_3,dist_source_pred_3,dist_source_actual_3,mag_source_3,magerr_source_3,source_density_1arcmin,dist_edge_left,dist_edge_right,dist_edge_bottom,dist_edge_top,preview_image_path,preview_image_file,status,detection_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)

                mcd.output_log_entry(path_logfile,'Inserting/updating multi-aperture photometry entry for detection_id {:d}...'.format(detection_id))
                query = "SELECT detection_id FROM detection_multiap_photometry WHERE detection_id={:d}".format(detection_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                row_detection = cursor.fetchone()
                if row_detection == None:
                    mcd.output_log_entry(path_logfile,'Inserting multi-aperture photometry entry for detection_id {:d}...'.format(detection_id))
                    query = "INSERT OR IGNORE INTO detection_multiap_photometry(detection_id,flux_t_05rKron,dflux_t_05rKron,mag_t_05rKron,magerr_t_05rKron,flux_t_10rKron,dflux_t_10rKron,mag_t_10rKron,magerr_t_10rKron,flux_t_20rKron,dflux_t_20rKron,mag_t_20rKron,magerr_t_20rKron,flux_t_30rKron,dflux_t_30rKron,mag_t_30rKron,magerr_t_30rKron,flux_t_40rKron,dflux_t_40rKron,mag_t_40rKron,magerr_t_40rKron,flux_t_50rKron,dflux_t_50rKron,mag_t_50rKron,magerr_t_50rKron) VALUES ({:d},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f})".format(detection_id,flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK,flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK,flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK,flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK,flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK,flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)
                else:
                    mcd.output_log_entry(path_logfile,'Updating multi-aperture photometry entry for detection_id {:d}...'.format(detection_id))
                    query = "UPDATE detection_multiap_photometry SET flux_t_05rKron={:.3f},dflux_t_05rKron={:.3f},mag_t_05rKron={:.3f},magerr_t_05rKron={:.3f},flux_t_10rKron={:.3f},dflux_t_10rKron={:.3f},mag_t_10rKron={:.3f},magerr_t_10rKron={:.3f},flux_t_20rKron={:.3f},dflux_t_20rKron={:.3f},mag_t_20rKron={:.3f},magerr_t_20rKron={:.3f},flux_t_30rKron={:.3f},dflux_t_30rKron={:.3f},mag_t_30rKron={:.3f},magerr_t_30rKron={:.3f},flux_t_40rKron={:.3f},dflux_t_40rKron={:.3f},mag_t_40rKron={:.3f},magerr_t_40rKron={:.3f},flux_t_50rKron={:.3f},dflux_t_50rKron={:.3f},mag_t_50rKron={:.3f},magerr_t_50rKron={:.3f} WHERE detection_id={:d}".format(flux_t_05rK,dflux_t_05rK,mag_t_05rK,magerr_t_05rK,flux_t_10rK,dflux_t_10rK,mag_t_10rK,magerr_t_10rK,flux_t_20rK,dflux_t_20rK,mag_t_20rK,magerr_t_20rK,flux_t_30rK,dflux_t_30rK,mag_t_30rK,magerr_t_30rK,flux_t_40rK,dflux_t_40rK,mag_t_40rK,magerr_t_40rK,flux_t_50rK,dflux_t_50rK,mag_t_50rK,magerr_t_50rK,detection_id)
                    mcd.output_log_entry(path_logfile,query)
                    cursor.execute(query)

                query = "UPDATE detections SET search_result_status='{:s}' WHERE detection_id={:d}".format(status,detection_id)
                mcd.output_log_entry(path_logfile,query)
                cursor.execute(query)
                                
            conn.commit() # Commit changes
            conn.close()  # Close connection to database file
    except Error as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: ingest_detection_data_file()'.format(detection_data_file_to_ingest))
        mcd.output_error_log_entry_nonstring(path_logfile,path_errorfile,e)
        keep_going = False
    return keep_going


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
        mcd.send_status_email('LN02_ingest_detection_data execution started','LN02_ingest_detection_data execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'LN02_ingest_detection_data')

        dir_detection_data = base_path + 'detection_data/'
        if not os.path.isdir(dir_detection_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Detection data directory {:s} not found'.format(dir_detection_data))
            mcd.send_status_email('LN02_ingest_detection_data execution failed','LN02_ingest_detection_data execution failed - Detection data directory {:s} not found.'.format(dir_detection_data))
            keep_going = False

    if keep_going:
        # Connect to database file
        os.chdir(dir_detection_data)
        for detection_data_file_to_ingest_gz in sorted(glob.glob('detection_data_*_toingest.txt.gz')):
            mcd.output_log_entry(path_logfile,'Ingesting {:s}...'.format(detection_data_file_to_ingest_gz))
            mcd.decompress_file_gzip(detection_data_file_to_ingest_gz)
            detection_data_file_to_ingest = detection_data_file_to_ingest_gz[:-3]
            keep_going = ingest_detection_data_file(detection_data_file_to_ingest,sqlite_file,keep_going,path_logfile,path_errorfile)
            if keep_going:
                ingested_filename = detection_data_file_to_ingest[:-13]+'_ingested.txt'
                os.rename(detection_data_file_to_ingest,ingested_filename)
                mcd.compress_file_gzip(ingested_filename)
            else:
                mcd.compress_file_gzip(detection_data_file_to_ingest)

    mcd.send_status_email('LN02_ingest_detection_data execution complete.','LN02_ingest_detection_data execution complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()


