import sys
import datetime
import glob, os, bz2, subprocess
import os.path
import macadamia_functions as mcd

def main():

    # Define filenames and paths
    if len(sys.argv)!=2:
        print('Usage:\n python3 pyt_MSC02a_stage_sdss_data_tocopy.py [thread_idx]\n')
    thread_idx = int(sys.argv[1])
    
    input_filename = '/data2/data_sdss_tocopy/file_list_{:02d}.txt'.format(thread_idx)
    destination_path = '/data2/data_sdss_tocopy/{:02d}/'.format(thread_idx)
    path_logfile = '/data2/data_sdss_tocopy/log_stage_sdss_data_tocopy_{:02d}_{:s}.txt'.format(thread_idx,datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    if not os.path.exists(destination_path):
        os.mkdir(destination_path)
    with open(input_filename,'r') as input_file, open(path_logfile,'w') as log_file:
        for line in input_file:
            source_files = '/data1/data_sdss_processed/' + line.strip()
            for file in sorted(glob.glob(source_files)):
                mcd.copy_file(file,destination_path)
                print('{:s} - Copying {:s} to {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),file,destination_path))
                log_file.write('{:s} - Copying {:s} to {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),file,destination_path))
    
    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
