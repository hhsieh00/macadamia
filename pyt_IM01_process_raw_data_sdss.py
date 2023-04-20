import sys
import datetime
import glob, os, bz2, subprocess
import os.path
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.table import Table
import sqlite3
from sqlite3 import Error
import smtplib, ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- BASIC DATA PROCESSING -- SDSS #####

def extract_extensions_sdss_data(mosaic_elem_table,base_path,dir_raw_data,dir_processed_data,sub_dir,data_dir,first_element,num_elements,keep_going,log_file):
    if keep_going:
        try:
            ### Set directory paths
            path_data = base_path+dir_raw_data+sub_dir+data_dir
            destination_dir = base_path+dir_processed_data+sub_dir+data_dir
            os.chdir(path_data)
            
            ### If directory is ready for processing, execute processing steps
            if os.path.isfile('00_download.complete'):
                if not os.path.isfile('00_extension_extraction.complete'):
                    # Check if extensions have already been extracted, and execute if not
                    print('{:s} - Processing data in {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),path_data))
                    log_file.write('{:s} - Processing data in {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),path_data))
                    # Create corresponding subdirectories in processed data directory if they do not already exist
                    log_file = mcd.create_directory(base_path+dir_processed_data,log_file)
                    log_file = mcd.create_directory(base_path+dir_processed_data+sub_dir,log_file)
                    log_file = mcd.create_directory(base_path+dir_processed_data+sub_dir+data_dir,log_file)
                    log_filename2 = path_data + '00_extension_extraction.inprogress'
                    data_table_filelinks_filename = path_data + 'data_table_filelinks.txt'
                    if os.path.isfile(log_filename2):
                        os.rename(log_filename2,path_data+'00_extension_extraction.inprogress.saved_{:s}'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S')))
                    with open(log_filename2,'w') as log_file2:
                        for data_filename in sorted(glob.glob('*.fz')):
                            if data_filename[6] != 'u':
                                if not os.path.isfile(destination_dir + data_filename[:-8] + '_000.fits.fz'):
                                    ### Decompress fz-compressed file, needed for reprocessing SDSS data
                                    print('{:s} - Processing {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                                    log_file2.write('{:s} - Processing {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                                    mcd.decompress_file_funpack(data_filename)
                                    ### Split multi-extension file into individual mosaic elements and save to processed directory
                                    keep_going,log_file2 = extract_extensions_sdss_file(mosaic_elem_table,data_filename[:-3],first_element,num_elements,path_data,destination_dir,keep_going,log_file2)
                                    mcd.compress_file_fpack(data_filename[:-3])
                                    keep_going = True
                        for data_filename in sorted(glob.glob('*.bz2')):
                            if data_filename[6] != 'u':
                                if not os.path.isfile(destination_dir + data_filename[:-9] + '_000.fits.fz'):
                                    ### Decompress bz2-compressed file, needed for SDSS data
                                    print('{:s} - Processing {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                                    log_file2.write('{:s} - Processing {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                                    mcd.decompress_file_bzip2(data_filename)
                                    ### Split multi-extension file into individual mosaic elements and save to processed directory
                                    keep_going,log_file2 = extract_extensions_sdss_file(mosaic_elem_table,data_filename[:-4],first_element,num_elements,path_data,destination_dir,keep_going,log_file2)
                                    mcd.compress_file_fpack(data_filename[:-4])
                                    keep_going = True
                    os.rename(log_filename2,path_data+'00_extension_extraction.complete')
                    keep_going,log_file = remove_intermediate_proc_files(path_data,keep_going,log_file)
                    print('{:s} - Extension extraction complete for {:s}.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                    log_file.write('{:s} - Extension extraction complete for {:s}.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename))
                    with open(destination_dir+'00_extension_extraction.complete','w') as f:
                        f.write('{:s} - Extension extraction complete.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            else:
                print('{:s} - Data downloading for {:s} not yet complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),path_data))
                log_file.write('{:s} - Data downloading for {:s} not yet complete.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),path_data))
        except:
            print('{:s} - Function failed: extract_extensions_sdss_data()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            print('{:s}'.format(base_path+dir_raw_data+sub_dir+data_dir))
            log_file.write('{:s} - Function failed: extract_extensions_sdss_data()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s}\n'.format(base_path+dir_raw_data+sub_dir+data_dir))
            mcd.send_status_email('IM01_process_raw_data_sdss failed.','{:s} - IM01_process_raw_data_sdss failed - extract_extensions_sdss_data().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        print('{:s} - Function skipped: extract_extensions_sdss_data()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Function skipped: extract_extensions_sdss_data()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return keep_going,log_file


def remove_intermediate_proc_files(path_data,keep_going,log_file):
    if keep_going:
        try:
            cmd_rm = 'rm'
            for inprogress_file in sorted(glob.glob('00_extension_extraction.inprogress*')):
                filepath = inprogress_file
                cmd = [cmd_rm,filepath]
                process = subprocess.call(cmd)
        except:
            print('{:s} - Function failed: remove_intermediate_proc_files()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - Function failed: remove_intermediate_proc_files()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM01_process_raw_data_sdss failed.','{:s} - IM01_process_raw_data_sdss failed - remove_intermediate_proc_files().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = False
    else:
        print('{:s} - Function skipped: remove_intermediate_proc_files()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Function skipped: remove_intermediate_proc_files()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

    return keep_going,log_file    
    
    
def extract_extensions_sdss_file(mosaic_elem_table,data_filename,first_element,num_elements,path_data,destination_dir,keep_going,log_file):
    # extract Extension 0 from multi-extension SDSS FITS files
    # return: writes individual extension files to working directory
    output_filename = ''
    if keep_going:
        try:
            print('{:s} - Splitting {:s} into {:d} individual mosaic elements...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename,num_elements))
            log_file.write('{:s} - Splitting {:s} into {:d} individual mosaic elements...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),data_filename,num_elements))
            multi_ext_filename = data_filename
            hdulist = fits.open(multi_ext_filename)
            #print(multi_ext_filename)
            #print(here)
            for idx in range(first_element,(first_element+num_elements)):
                output_filename = multi_ext_filename[:-5] + '_{:03d}.fits'.format(idx)
                hdr1 = getheader(multi_ext_filename,0)
                radecsys = hdr1['RADECSYS']
                hdr1['RADESYSA'] = radecsys
                del hdr1['RADECSYS']
                hdr1['INSTRUME'] = 'SDSS'
                extension_data = hdulist[idx].data
                if os.path.isfile(output_filename):
                    print('{:s} - Output fits file {:s} already exists...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),output_filename))
                    log_file.write('{:s} - Output fits file {:s} already exists...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),output_filename))
                else:
                    print('{:s} - Writing fits file {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),output_filename))
                    log_file.write('{:s} - Writing fits file {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),output_filename))
                    fits.writeto(output_filename,extension_data,hdr1)
            hdulist.close()
            mcd.compress_file_fpack(output_filename)
            os.rename(path_data+output_filename+'.fz',destination_dir+output_filename+'.fz')
            print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            keep_going = True
        except Error as e:
            print('{:s} - Function failed: extract_extensions_sdss_file()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('\n{:s} - Function failed: extract_extensions_sdss_file()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM01_process_raw_data_sdss failed.','{:s} - IM01_process_raw_data_sdss failed - extract_extensions_sdss_file().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            print(e)
            keep_going = False
    else:
        print('{:s} - Function skipped: extract_extensions_sdss_file()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Function skipped: extract_extensions_sdss_file()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return keep_going,log_file


##### FUNCTION DEFINITIONS -- DATABASE -- RETRIEVE DATA #####

# Retrieve instrument_id corresponding to a given instrument name
def retrieve_instrument_id(sqlite_file,instrument_name,keep_going,log_file):
    instrument_id = int(0)
    if keep_going:
        print('{:s} - Retrieving instrument ID for {:s} from database...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
        log_file.write('{:s} - Retrieving instrument ID for {:s} from database...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            cursor.execute("SELECT instrument_id FROM instruments WHERE instrument_name='{:s}'".format(instrument_name))
            log_file.write("{:s} - QUERY: SELECT instrument_id FROM instruments WHERE instrument_name='{:s}'\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
            row = cursor.fetchone()
            if row != None:
                instrument_id = int(row[0])
                print('{:s} - instrument_id={:d}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_id))
                log_file.write('{:s} - instrument_id={:d}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_id))
                keep_going = True
            else:
                print("{:s} - Instrument '{:s}' not found.".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                log_file.write("{:s} - Instrument '{:s}' not found.\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Error as e:
            print('{:s} - Function failed: retrieve_instrument_id()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - Function failed: retrieve_instrument_id()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM01_process_raw_data_sdss failed.','{:s} - IM01_process_raw_data_sdss failed - retrieve_instrument_id().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            print(e)
            keep_going = False
    else:
        print('{:s} - Function skipped: retrieve_instrument_id()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Function skipped: retrieve_instrument_id()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return(instrument_id,keep_going,log_file)


# Retrieve mosaic_element_id corresponding to a given mosaic element for a given instrument
def retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,keep_going,log_file):
    mosaic_element_id = -1
    if keep_going:
        print('{:s} - Retrieving mosaic element ID from database...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Retrieving mosaic element ID from database...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            cursor.execute("SELECT mosaic_element_id FROM mosaic_elements WHERE instrument_id={:d} AND mosaic_element_num={:d}".format(instrument_id,mosaic_element_num))
            log_file.write("{:s} - SELECT mosaic_element_id FROM mosaic_elements WHERE instrument_id={:d} AND mosaic_element_num={:d}\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),instrument_id,mosaic_element_num))
            row = cursor.fetchone()
            if row != None:
                mosaic_element_id = int(row[0])
                print("{:s} - Mosaic element {:d} (mosaic_element_id={:d}) for {:s} instrument found.".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,mosaic_element_id,instrument_name))
                log_file.write("{:s} - Mosaic element {:d} (mosaic_element_id={:d}) for {:s} instrument found.\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,mosaic_element_id,instrument_name))
                print('{:s} - mosaic_element_id={:d}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_id))
                log_file.write('{:s} - mosaic_element_id={:d}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_id))
                keep_going = True
            else:
                print("{:s} - Mosaic element {:d} for {:s} instrument not found.".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                log_file.write("{:s} - Mosaic element {:d} for {:s} instrument not found.\n".format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),mosaic_element_num,instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Error as e:
            print('{:s} - Function failed: retrieve_mosaic_element_id()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            log_file.write('{:s} - Function failed: retrieve_mosaic_element_id()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            mcd.send_status_email('IM01_process_raw_data_sdss failed.','{:s} - IM01_process_raw_data_sdss failed - retrieve_mosaic_element_id().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            print(e)
            keep_going = False
    else:
        print('{:s} - Function skipped: retrieve_mosaic_element_id()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Function skipped: retrieve_mosaic_element_id()\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return(mosaic_element_id,keep_going,log_file)

    
def main():

    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_IM01_process_raw_data_sdss.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = sys.argv[2]

    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)
    
    if keep_going:
        mcd.send_status_email('IM01_process_raw_data_sdss execution started','{:s} - IM01_process_raw_data_sdss execution started.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        
        # Create and initialize log file
        path_logfile = mcd.initialize_log_file(base_path,'IM01_process_raw_data_sdss')

        os.chdir(base_path)
        dir_raw_data       = 'data_sdss_raw/'
        dir_processed_data = 'data_sdss_processed/'
        with open(path_logfile,'a') as log_file:
            if not os.path.isdir(dir_raw_data):
                print('{:s} - Raw data directory {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),dir_raw_data))
                log_file.write('{:s} - Raw data directory {:s} not found.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),dir_raw_data))
                mcd.send_status_email('IM01_process_raw_data_sdss execution failed','{:s} - IM01_process_raw_data_sdss execution failed - Raw data directory {:s} not found.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),dir_raw_data))
                keep_going = False
    
    if keep_going:
        with open(path_logfile,'a') as log_file:
            ### SDSS-specific data ###
            instrument_name = 'SDSS Imaging Camera'
            first_element = 0
            num_elements  = 1
            mosaic_elem_datarows = [(0,-1)]  # element index number and dummy value for element_id

            # Retrieve telescope/instrument data
            instrument_id,keep_going,log_file = retrieve_instrument_id(sqlite_file,instrument_name,keep_going,log_file)
            mosaic_elem_table    = Table(rows=mosaic_elem_datarows,names=('elem_idx','element_id'))
            for idx in range(0,num_elements):
                mosaic_element_num = mosaic_elem_table[idx]['elem_idx']
                mosaic_element_id,keep_going,log_file = retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,keep_going,log_file)
                mosaic_elem_table[idx]['element_id'] = mosaic_element_id
    
            # Process data
            mcd.create_directory(dir_processed_data,log_file)
            os.chdir(base_path+dir_raw_data)
            for sub_dir in sorted(glob.glob('*/')):
                os.chdir(base_path+dir_raw_data+sub_dir)
                for data_dir in sorted(glob.glob('*/')):
                    keep_going = True
                    # Extract data extensions for all data in directory
                    keep_going,log_file = extract_extensions_sdss_data(mosaic_elem_table,base_path,dir_raw_data,dir_processed_data,sub_dir,data_dir,first_element,num_elements,keep_going,log_file)

        mcd.send_status_email('IM01_process_raw_data_sdss execution complete','{:s} - IM01_process_raw_data_sdss execution complete.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
