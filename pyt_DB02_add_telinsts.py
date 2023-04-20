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
    sqlite_file = base_path + sys.argv[2]
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('DB02_add_telinsts execution started','{:s} - DB02_add_telinsts execution started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
        # Create and initialize log file
        path_logfile = mcd.initialize_log_file(base_path,'DB02_add_telinsts')
        
        with open(path_logfile,'a') as log_file:
            telescope_file   = base_path + 'telinst_data/list_telescopes.dat'
            instrument_file  = base_path + 'telinst_data/list_instruments.dat'
            mosaicelems_file = base_path + 'telinst_data/list_mosaic_elements.dat'
        
            if not os.path.isfile(telescope_file):
                print('{:s} - Telescope data file {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_file))
                log_file.write('{:s} - Telescope data file {:s} not found.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_file))
                mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Telescope data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),telescope_file))
                keep_going = False
            if not os.path.isfile(instrument_file):
                print('{:s} - Instrument data file {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_file))
                log_file.write('{:s} - Instrument data file {:s} not found.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_file))
                mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Instrument data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),instrument_file))
                keep_going = False
            if not os.path.isfile(mosaicelems_file):
                print('{:s} - Mosaic elements data file {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaicelems_file))
                log_file.write('{:s} - Mosaic elements data file {:s} not found.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaicelems_file))
                mcd.send_status_email('DB02_add_telinsts execution failed','{:s} - DB02_add_telinsts execution failed - Mosaic elements data file {:s} not found.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),mosaicelems_file))
                keep_going = False
    
    if keep_going:
        with open(path_logfile,'a') as log_file:
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
                        cursor.execute("SELECT telescope_id FROM telescopes WHERE telescope_name='{:s}'".format(telescope_data[0]))
                        row = cursor.fetchone()
                        if row == None: # If no telescope entry exists, add new entry
                            cursor.execute("INSERT INTO telescopes(telescope_name,observatory_code,aperture_size,latitude,longitude,elevation) VALUES ('{:s}','{:s}',{:5.2f},{:9.5f},{:9.5f},{:7.1f})".format(telescope_data[0],telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5])))
                            log_file.write("{:s} - QUERY: INSERT INTO telescopes(telescope_name,observatory_code,aperture_size,latitude,longitude,elevation) VALUES ('{:s}','{:s}',{:5.2f},{:9.5f},{:9.5f},{:7.1f})\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[0],telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5])))
                            print('{:s} - Telescope data added for {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[0]))
                            log_file.write('{:s} - Telescope data added for {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[0]))
                        else:
                            cursor.execute("UPDATE telescopes SET observatory_code='{:s}',aperture_size={:5.2f},latitude={:9.5f},longitude={:9.5f},elevation={:7.1f} WHERE telescope_name='{:s}'".format(telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5]),telescope_data[0]))
                            log_file.write("{:s} - UPDATE telescopes SET observatory_code='{:s}',aperture_size={:5.2f},latitude={:9.5f},longitude={:9.5f},elevation={:7.1f} WHERE telescope_name='{:s}'\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[1],float(telescope_data[2]),float(telescope_data[3]),float(telescope_data[4]),float(telescope_data[5]),telescope_data[0]))
                            print('{:s} - Telescope data updated for {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[0]))
                            log_file.write('{:s} - Telescope data updated for {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_data[0]))
                        line_count += 1
            print('{:s} - >>> Telescope data added/updated.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Telescope data added/updated.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

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
                        cursor.execute("SELECT telescope_id FROM telescopes WHERE telescope_name='{:s}'".format(telescope_name))
                        row = cursor.fetchone()
                        if row == None:
                            print("{:s} - Telescope '{:s}' for instrument '{:s}' not found!".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_name,instrument_name))
                            log_file.write("{:s} - Telescope '{:s}' for instrument '{:s}' not found!\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_name,instrument_name))
                        else:
                            telescope_id = row[0]
                            # If telescope is found, check if instrument entry exists
                            cursor.execute("SELECT instrument_id FROM instruments WHERE instrument_name='{:s}' and telescope_id={:d}".format(instrument_name,telescope_id))
                            row = cursor.fetchone()
                            if row == None: # If no instrument entry exists, add new entry
                                cursor.execute("INSERT INTO instruments(telescope_id,instrument_name) VALUES ({:d},'{:s}')".format(telescope_id,instrument_name))
                                log_file.write("{:s} - QUERY: INSERT INTO instruments(telescope_id,instrument_name) VALUES ({:d},'{:s}')\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_id,instrument_name))
                                print('{:s} - Instrument data added for {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                                log_file.write('{:s} - Instrument data added for {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                            else: # If instrument entry exists, update entry
                                instrument_id = row[0]
                                cursor.execute("UPDATE instruments SET telescope_id={:d},instrument_name='{:s}' WHERE instrument_id={:d}".format(telescope_id,instrument_name,instrument_id))
                                log_file.write("{:s} - QUERY: UPDATE instruments SET telescope_id={:d},instrument_name='{:s}' WHERE instrument_id={:d}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),telescope_id,instrument_name,instrument_id))
                                print('{:s} - Instrument data updated for {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                                log_file.write('{:s} - Instrument data updated for {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                        line_count += 1
            print('{:s} - >>> Instrument data added/updated.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Instrument data added/updated.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

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
                        cursor.execute("SELECT inst.instrument_id FROM instruments AS inst INNER JOIN telescopes AS tel ON inst.telescope_id=tel.telescope_id WHERE inst.instrument_name='{:s}' and tel.telescope_name='{:s}'".format(instrument_name,telescope_name))
                        row = cursor.fetchone()
                        if row == None:
                            print("{:s} - Instrument '{:s}' on telescope '{:s}' not found!".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name,telescope_name))
                            log_file.write("{:s} - Instrument '{:s}' on telescope '{:s}' not found!\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name,telescope_name))
                        else:
                            instrument_id = row[0]
                            # If telescope and instrument are found, check if mosaic entry exists
                            cursor.execute("SELECT mosaic_element_id FROM mosaic_elements WHERE mosaic_element_num={:d} AND instrument_id={:d}".format(mosaic_element_num,instrument_id))
                            row = cursor.fetchone()
                            if row == None: # If no mosaic element entry exists, add new entry
                                cursor.execute("INSERT INTO mosaic_elements(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y) VALUES ({:d},{:d},{:f},{:f},{:d},{:d},{:f},{:f})".format(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y))
                                log_file.write("{:s} - QUERY: INSERT INTO mosaic_elements(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y) VALUES ({:d},{:d},{:f},{:f},{:d},{:d},{:f},{:f})\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y))
                                print('{:s} - Mosaic element data added for Element {:d} of the {:s} instrument...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                                log_file.write('{:s} - Mosaic element data added for Element {:d} of the {:s} instrument...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                            else: # If mosaic element entry exists, update entry
                                mosaic_element_id = row[0]
                                cursor.execute("UPDATE mosaic_elements SET instrument_id='{:d}',mosaic_element_num='{:d}',gain={:f},read_noise={:f},npix_x={:d},npix_y={:d},pixel_scale_x={:f},pixel_scale_y={:f} WHERE mosaic_element_id={:d}".format(instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y,mosaic_element_id))
                                log_file.write("{:s} - QUERY: UPDATE mosaic_elements SET instrument_id='{:d}',mosaic_element_num='{:d}',gain={:f},read_noise={:f},npix_x={:d},npix_y={:d},pixel_scale_x={:f},pixel_scale_y={:f} WHERE mosaic_element_id={:d}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_id,mosaic_element_num,gain,read_noise,npix_x,npix_y,pixel_scale_x,pixel_scale_y,mosaic_element_id))
                                print('{:s} - Mosaic element data updated for Element {:d} of the {:s} instrument...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                                log_file.write('{:s} - Mosaic element data updated for Element {:d} of the {:s} instrument...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                        line_count += 1
                    
            print('{:s} - >>> Mosaic element data added/updated.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Mosaic element data added/updated.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
            # Committing changes and closing the connection to the database file
            conn.commit()
            conn.close()
            print('{:s} - >>> Database closed.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - >>> Database closed.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

        mcd.send_status_email('DB02_add_telinsts execution complete','{:s} - DB02_add_telinsts execution complete.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()

