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

### Edit History ###
# 2021-06-28: check extraction completion statuses of individual files rather than entire directories

##### FUNCTION DEFINITIONS -- BASIC DATA PROCESSING -- MEGACAM #####

def extract_extensions_megacam_data(mosaic_elem_table,base_path,dir_raw_data,dir_processed_data,data_dir,first_element,num_elements,thread_idx,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            ### Set directory paths
            path_data = base_path+dir_raw_data+data_dir
            destination_dir = base_path+dir_processed_data+data_dir
            os.chdir(path_data)
            
            # Check if extensions have already been extracted, and execute if not
            mcd.output_log_entry(path_logfile,'Processing data in {:s}...'.format(path_data))
            # Create corresponding subdirectories in processed data directory if they do not already exist
            mcd.create_directory(base_path+dir_processed_data,path_logfile,path_errorfile)
            mcd.create_directory(base_path+dir_processed_data+data_dir,path_logfile,path_errorfile)
            for data_filename in sorted(glob.glob('*.fz')):
                # check that last extracted element of current raw file does not exist in processed data directory (meaning file still needs processing)
                if not os.path.isfile(destination_dir+data_filename[:-8]+'_{:03d}.fits.fz'.format(first_element+num_elements-1)):
                    mcd.output_log_entry(path_logfile,'Processing {:s}...'.format(data_filename))
                    # if current file is partially extracted, delete extracted elements and re-process
                    if os.path.isfile(destination_dir+data_filename[:-8]+'_{:03d}.fits.fz'.format(first_element)):
                        mcd.output_log_entry(path_logfile,'{:s} partially processed; deleting old files and re-processing...'.format(data_filename))
                        for idx in range(first_element,first_element+num_elements-1):
                            if os.path.isfile(destination_dir+data_filename[:-8]+'_{:03d}.fits.fz'.format(idx)):
                                os.remove(destination_dir+data_filename[:-8]+'_{:03d}.fits.fz'.format(idx))
                    ### Decompress fz-compressed file
                    mcd.decompress_file_funpack(data_filename)
                    ### Split multi-extension file into individual mosaic elements and save to processed directory
                    keep_going = extract_extensions_megacam_file(mosaic_elem_table,data_filename[:-3],first_element,num_elements,path_data,destination_dir,thread_idx,keep_going,path_logfile,path_errorfile)
                    mcd.compress_file_fpack(data_filename[:-3])
                    keep_going = True
            keep_going = remove_intermediate_proc_files(path_data,thread_idx,keep_going,path_logfile,path_errorfile)
            mcd.output_log_entry(path_logfile,'Extension extraction complete for {:s}.'.format(data_filename))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: extract_extensions_megacam_data()'.format(base_path+dir_raw_data+data_dir))
            print(e)
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} failed.'.format(thread_idx),'IM01_process_raw_data_megacam_{:02d} failed - extract_extensions_megacam_data().'.format(thread_idx))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: extract_extensions_megacam_data()')
    return keep_going


def remove_intermediate_proc_files(path_data,thread_idx,keep_going,path_logfile,path_errorfile):
    if keep_going:
        try:
            for inprogress_file in sorted(glob.glob('00_extension_extraction.inprogress*')):
                filepath = inprogress_file
                os.remove(filepath)
        except:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: remove_intermediate_proc_files()'.format(path_data))
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} failed.'.format(thread_idx),'{:s} - IM01_process_raw_data_megacam_{:02d} failed - remove_intermediate_proc_files().'.format(path_data,thread_idx))
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: remove_intermediate_proc_files()')
    return keep_going
    

def extract_extensions_megacam_file(mosaic_elem_table,data_filename,first_element,num_elements,path_data,destination_dir,thread_idx,keep_going,path_logfile,path_errorfile):
    # extract Extensions 1-40 from multi-extension CFHT MegaCam files
    # return: writes individual extension files to working directory
    output_filename = ''
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Splitting {:s} into {:d} individual mosaic elements...'.format(data_filename,num_elements))
            multi_ext_filename = data_filename
            hdulist = fits.open(multi_ext_filename)
            num_elements = len(hdulist) - 1
            for idx in range(first_element,(first_element+num_elements)):
                output_filename = multi_ext_filename[:-5] + '_{:03d}.fits'.format(idx)
                #hdr1 = getheader(multi_ext_filename,0)
                #radecsys = hdr1['RADECSYS']
                #hdr1['RADESYSA'] = radecsys
                #del hdr1['RADECSYS']
                extension_data = hdulist[idx].data
                hdr1 = hdulist[idx].header
                hdr1['INSTRUME'] = 'MegaCam'
                if os.path.isfile(output_filename) or os.path.isfile(output_filename+'.fz'):
                    mcd.output_log_entry(path_logfile,'Output fits file {:s} already exists...'.format(output_filename))
                else:
                    mcd.output_log_entry(path_logfile,'Writing fits file {:s}...'.format(output_filename))
                    fits.writeto(output_filename,extension_data,hdr1)
                    mcd.compress_file_fpack(output_filename)
                    os.rename(path_data+output_filename+'.fz',destination_dir+output_filename+'.fz')
            hdulist.close()
            mcd.output_log_entry(path_logfile,'Splitting {:s} into individual mosaic elements done.'.format(data_filename))
            keep_going = True
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: extract_extensions_megacam_file()'.format(path_data+data_filename))
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} failed.'.format(thread_idx),'{:s} - IM01_process_raw_data_megacam_{:02d} failed - extract_extensions_megacam_file().'.format(path_data+data_filename,thread_idx))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: extract_extensions_megacam_file()')
    return keep_going


##### FUNCTION DEFINITIONS -- DATABASE -- RETRIEVE DATA #####

def retrieve_instrument_id(sqlite_file,instrument_name,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve instrument_id corresponding to a given instrument name
    instrument_id = int(0)
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
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Instrument {:s} not found.'.format(instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} - Function failed: retrieve_instrument_id()'.format(instrument_name))
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} failed.'.format(thread_idx),'{:s} - IM01_process_raw_data_megacam_{:02d} failed - retrieve_instrument_id().'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),thread_idx))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_instrument_id()')
    return(instrument_id,keep_going)


def retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile):
    # Retrieve mosaic_element_id corresponding to a given mosaic element for a given instrument
    mosaic_element_id = -1
    if keep_going:
        mcd.output_log_entry(path_logfile,'Retrieving mosaic element ID ({:s}, {:d}) from database...'.format(instrument_name,mosaic_element_num))
        try:
            conn = mcd.create_connection(sqlite_file)  # Open connection to database file
            cursor = conn.cursor()
            query = 'SELECT mosaic_element_id FROM mosaic_elements WHERE instrument_id={:d} AND mosaic_element_num={:d}'.format(instrument_id,mosaic_element_num)
            cursor.execute(query)
            mcd.output_log_entry(path_logfile,query)
            row = cursor.fetchone()
            if row != None:
                mosaic_element_id = int(row[0])
                mcd.output_log_entry(path_logfile,'Mosaic element {:d} (mosaic_element_id={:d}) for {:s} instrument found.'.format(mosaic_element_num,mosaic_element_id,instrument_name))
                keep_going = True
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'Mosaic element {:d} for {:s} instrument not found.'.format(mosaic_element_num,instrument_name))
                keep_going = False
            conn.close()  # Close connection to database file
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'{:s} Ext {:d} - Function failed: retrieve_mosaic_element_id()'.format(instrument_name,mosaic_element_num))
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} failed.'.format(thread_idx),'{:s} Ext {:d} - IM01_process_raw_data_megacam_{:02d} failed - retrieve_mosaic_element_id().'.format(instrument_name,mosaic_element_num,thread_idx))
            print(e)
            keep_going = False
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: retrieve_mosaic_element_id()')
    return mosaic_element_id,keep_going

    
def main():

    # Define filenames and paths
    if len(sys.argv)!=4:
        print('Usage:\n python3 pyt_IM01_process_raw_data_megacam.py [base_path] [sqlite_file] [thread_idx]')
        print(" (Trailing '/' needed in path specification)")
        print(" (Full path needed for SQLite file)")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    thread_idx  = int(sys.argv[3])

    # Validate input parameters
    keep_going = True
    if base_path[-1:] != '/': base_path = base_path + '/'
    if not os.path.isdir(base_path):
        print('Directory {:s} not found.'.format(base_path))
        keep_going = False
    
    if keep_going:
        mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} execution started'.format(thread_idx),'IM01_process_raw_data_megacam_{:02d} execution started.'.format(thread_idx))
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'IM01_process_raw_data_megacam_{:02d}'.format(thread_idx))

        os.chdir(base_path)
        dir_raw_data       = 'data_megacam_raw/'
        dir_processed_data = 'data_megacam_processed/'
        if not os.path.isdir(dir_raw_data):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Raw data directory {:s} not found.')
            mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} execution failed'.format(thread_idx),'IM01_process_raw_data_megacam_{:02d} execution failed - Raw data directory {:s} not found.'.format(thread_idx,dir_raw_data))
            keep_going = False
    
    if keep_going:
        os.chdir(base_path)
    
        ### MegaCam-specific data ###
        instrument_name = 'MegaCam'
        first_element = 1
        num_elements  = 40
        mosaic_elem_datarows = [[0 for idx2 in range(2)] for idx1 in range(num_elements)]
        #mosaic_elem_datarows = [(0,-1)]  # element index number and dummy value for element_id

        # Retrieve telescope/instrument data
        instrument_id,keep_going = retrieve_instrument_id(sqlite_file,instrument_name,thread_idx,keep_going,path_logfile,path_errorfile)
        mosaic_elem_table = Table(rows=mosaic_elem_datarows,names=('elem_idx','element_id'))
        for idx in range(0,num_elements):
            mosaic_elem_table[idx]['elem_idx'] = first_element+idx
            mosaic_element_num = mosaic_elem_table[idx]['elem_idx']
            mosaic_element_id,keep_going = retrieve_mosaic_element_id(sqlite_file,instrument_id,instrument_name,mosaic_element_num,thread_idx,keep_going,path_logfile,path_errorfile)
            mosaic_elem_table[idx]['element_id'] = mosaic_element_id
    
        # Process data
        mcd.create_directory(dir_processed_data,path_logfile,path_errorfile)
        os.chdir(base_path+dir_raw_data)
        for data_dir in sorted(glob.glob('ut*{:d}/'.format(thread_idx))):
            os.chdir(base_path+dir_raw_data+data_dir)
            keep_going = True
            # Extract data extensions for all data in directory
            keep_going = extract_extensions_megacam_data(mosaic_elem_table,base_path,dir_raw_data,dir_processed_data,data_dir,first_element,num_elements,thread_idx,keep_going,path_logfile,path_errorfile)

        mcd.send_status_email('IM01_process_raw_data_megacam_{:02d} execution complete'.format(thread_idx),'IM01_process_raw_data_megacam_{:02d} execution complete.'.format(thread_idx))

    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()
