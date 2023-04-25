import re
import sys
import time
import datetime
import glob, os, bz2, subprocess
import csv
import os.path
import sqlite3
from sqlite3 import Error
import macadamia_functions as mcd


def main():
    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_add_telinsts.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DB02_add_telinsts execution started','{:s} - DB02_add_telinsts execution started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
        # Create and open log file
        dir_logfiles = base_path + 'log_files/'
        if not os.path.isdir(dir_logfiles): os.mkdir(dir_logfiles)
        path_logfile = dir_logfiles + 'log_{:s}_DB02_add_telinsts.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
        with open(path_logfile,'w') as log_file:
            log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
        telescope_file   = base_path + 'telinst_data/list_telescopes.dat'
        instrument_file  = base_path + 'telinst_data/list_instruments.dat'
        mosaicelems_file = base_path + 'telinst_data/list_mosaic_elements.dat'
        
        if not os.path.isfile(telescope_file):
            mcd.output_log_entry(path_logfile,'Telescope data file {:s} not found.'.format(telescope_file))
            mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Telescope data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),telescope_file))
            keep_going = False
        if not os.path.isfile(instrument_file):
            mcd.output_log_entry(path_logfile,'Instrument data file {:s} not found.'.format(instrument_file))
            mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Instrument data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),instrument_file))
            keep_going = False
        if not os.path.isfile(mosaicelems_file):
            mcd.output_log_entry(path_logfile,'Mosaic elements data file {:s} not found.'.format(mosaicelems_file))
            mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Mosaic elements data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),mosaicelems_file))
            keep_going = False
    
    if keep_going:
        # Connect to database file
        conn = mcd.create_connection(sqlite_file)
        cursor = conn.cursor()

        with open(telescope_file) as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=',')
            line_count = 0
            for telescope_data in csv_reader:
                if line_count == 0:
                    line_count += 1
                else:
                    query = "SELECT telescope_id FROM telescopes WHERE telescope_name='{:s}'".format(telescope_data[0])
                    cursor.execute(query)
                    mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                    row = cursor.fetchone()
                    if row == None: # If no telescope entry exists, add new entry
                        query = "INSERT INTO telescopes(telescope_name,observatory_code,aperture_size,latitude,longitude,elevation) VALUES ('{:s}','{:s}',{:5.2f},{:9.5f},{:9.5f},{:7.1f})".format(telescope_data[0],telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5]))
                        cursor.execute(query)
                        mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                        mcd.output_log_entry(path_logfile,'Telescope data added for {:s}...'.format(telescope_data[0]))
                    else:
                        query = "UPDATE telescopes SET observatory_code='{:s}',aperture_size={:5.2f},latitude={:9.5f},longitude={:9.5f},elevation={:7.1f} WHERE telescope_name='{:s}'".format(telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5]),telescope_data[0])
                        cursor.execute(query)
                        mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                        mcd.output_log_entry(path_logfile,'Telescope data updated for {:s}...'.format(telescope_data[0]))
                    line_count += 1
        mcd.output_log_entry(path_logfile,'>>> Telescope data added/updated.')

        with open(instrument_file) as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=',')
            line_count = 0
            for instrument_data in csv_reader:
                if line_count == 0:  # Skip first header line
                    line_count += 1
                else:
                    instrument_name = instrument_data[0]
                    telescope_name  = instrument_data[1]
                    # Check if corresponding telescope listed for instrument is in telescope table
                    query = "SELECT telescope_id FROM telescopes WHERE telescope_name='{:s}'".format(telescope_name)
                    cursor.execute(query)
                    mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                    row = cursor.fetchone()
                    if row == None:
                        mcd.output_log_entry(path_logfile,"Telescope '{:s}' for instrument '{:s}' not found.".format(telescope_name,instrument_name))
                    else:
                        telescope_id = row[0]
                        # If telescope is found, check if instrument entry exists
                        query = "SELECT instrument_id FROM instruments WHERE instrument_name='{:s}' and telescope_id={:d}".format(instrument_name,telescope_id)
                        cursor.execute(query)
                        mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                        row = cursor.fetchone()
                        if row == None: # If no instrument entry exists, add new entry
                            query = "INSERT INTO instruments(telescope_id,instrument_name) VALUES ({:d},'{:s}')".format(telescope_id,instrument_name)
                            cursor.execute(query)
                            mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                            mcd.output_log_entry(path_logfile,'Instrument data added for {:s}...'.format(instrument_name))
                        else: # If instrument entry exists, update entry
                            instrument_id = row[0]
                            query = "UPDATE instruments SET telescope_id={:d},instrument_name='{:s}' WHERE instrument_id={:d}".format(telescope_id,instrument_name,instrument_id)
                            cursor.execute(query)
                            mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                            mcd.output_log_entry(path_logfile,'Instrument data updated for {:s}...'.format(instrument_name))
                    line_count += 1
        mcd.output_log_entry(path_logfile,'>>> Instrument data added/updated.')

        with open(mosaicelems_file) as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=',')
            line_count = 0
            for mosaic_element_data in csv_reader:
                if line_count == 0:  # Skip first header line
                    line_count += 1
                else:
                    mosaic_element_num = int(mosaic_element_data[0])
                    instrument_name    = mosaic_element_data[1]
                    telescope_name     = mosaic_element_data[2]
                    gain               = float(mosaic_element_data[3])
                    read_noise         = float(mosaic_element_data[4])
                    npix_x             = int(mosaic_element_data[5])
                    npix_y             = int(mosaic_element_data[6])
                    pixel_scale_x      = float(mosaic_element_data[7])
                    pixel_scale_y      = float(mosaic_element_data[8])
                    # Check if corresponding telescope and instrument listed for mosaic element are in corresponding tables
                    query = "SELECT inst.instrument_id FROM instruments AS inst INNER JOIN telescopes AS tel ON inst.telescope_id=tel.telescope_id WHERE inst.instrument_name='{:s}' and tel.telescope_name='{:s}'".format(instrument_name,telescope_name)
                    cursor.execute(query)
                    mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                    row = cursor.fetchone()
                    if row == None:
                        mcd.output_log_entry(path_logfile,"Instrument '{:s}' on telescope '{:s}' not found.".format(instrument_name,telescope_name))
                    else:
                        instrument_id = row[0]
                        # If telescope and instrument are found, check if mosaic entry exists
                        query = "SELECT mosaic_element_id FROM mosaic_elements WHERE mosaic_element_num={:d} AND instrument_id={:d}".format(mosaic_element_num,instrument_id)
                        cursor.execute(query)
                        mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                        row = cursor.fetchone()
                        if row == None: # If no mosaic element entry exists, add new entry
                            query = "INSERT INTO mosaic_elements(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y) VALUES ({:d},{:d},{:f},{:f},{:d},{:d},{:f},{:f})".format(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y)
                            cursor.execute(query)
                            mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                            mcd.output_log_entry(path_logfile,"Mosaic element data added for Element {:d} of the {:s} instrument...".format(mosaic_element_num,instrument_name))
                        else: # If mosaic element entry exists, update entry
                            mosaic_element_id = row[0]
                            query = "UPDATE mosaic_elements SET instrument_id='{:d}',mosaic_element_num='{:d}',gain={:f},read_noise={:f},npix_x={:d},npix_y={:d},pixel_scale_x={:f},pixel_scale_y={:f} WHERE mosaic_element_id={:d}".format(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y,mosaic_element_id)
                            cursor.execute(query)
                            mcd.output_log_entry(path_logfile,'QUERY: {:s}'.format(query))
                            mcd.output_log_entry(path_logfile,"Mosaic element data updated for Element {:d} of the {:s} instrument...".format(mosaic_element_num,instrument_name))
                    line_count += 1
                    
        mcd.output_log_entry(path_logfile,">>> Mosaic element data added/updated.")
        
        # Committing changes and closing the connection to the database file
        conn.commit()
        conn.close()
        mcd.output_log_entry(path_logfile,">>> Database closed.")

        mcd.send_status_email('DB02_add_telinsts execution complete','{:s} - DB02_add_telinsts execution complete.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    return None


if __name__ == '__main__':
    main()

