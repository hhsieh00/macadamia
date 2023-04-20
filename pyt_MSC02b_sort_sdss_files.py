import sys
import datetime
import glob, os, bz2, subprocess
import os.path
import macadamia_functions as mcd

def sort_sdss_file(fits_filename,log_file):
    data_dir         = int(fits_filename[8:14])
    data_subdir      = int(fits_filename[15:16])
    data_dir_path    = '/home/hhsieh/macadamia/data_sdss_processed/{:d}/'.format(data_dir)
    data_subdir_path = '/home/hhsieh/macadamia/data_sdss_processed/{:d}/{:d}/'.format(data_dir,data_subdir)
    if not os.path.exists(data_dir_path):
        os.mkdir(data_dir_path)
    if not os.path.exists(data_subdir_path):
        os.mkdir(data_subdir_path)
    fits_filestem = fits_filename[:21]
    print('{:s} - Moving {:s}* to {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),fits_filestem,data_subdir_path))
    log_file.write('{:s} - Moving {:s}* to {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),fits_filestem,data_subdir_path))
    for file in sorted(glob.glob(fits_filestem+'*')):
        mcd.move_file(file,data_subdir_path)
    return log_file


def main():

    # Define filenames and paths
    if len(sys.argv)!=2:
        print('Usage:\n python3 pyt_MSC02b_sort_sdss_files.py [staged_data_directory]\n')
    staged_data_directory = sys.argv[1]
    
    path_logfile = '/home/hhsieh/macadamia/log_sort_sdss_files_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    os.chdir(staged_data_directory)
    with open(path_logfile,'w') as log_file:
        for fits_filename in sorted(glob.glob('*.fits.fz')):
            log_file = sort_sdss_file(fits_filename,log_file)
            
    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
