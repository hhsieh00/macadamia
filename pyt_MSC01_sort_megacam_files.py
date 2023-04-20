import sys
import datetime
import glob, os, bz2, subprocess
import os.path
from astropy.io import fits
from astropy.io.fits import getheader
import smtplib, ssl
import macadamia_functions as mcd

def extract_date_megacam(fits_filename,keep_going):
    date_obs = ''
    hdulist = fits.open(fits_filename)
    hdr  = hdulist[0].header
    date_obs = hdr['DATE-OBS']
    return date_obs,keep_going

def main():

    # Define filenames and paths
    if len(sys.argv)!=2:
        print('Usage:\n python3 pyt_MSC01_sort_megacam_files.py [base_path]\n')
    base_path = sys.argv[1]
    
    keep_going = True
    
    # Validate input parameters
    if not os.path.isdir(base_path):
        print('Directory {:s} not found.'.format(base_path))
        keep_going = False
    
    if keep_going:
        mcd.send_status_email('MSC01_sort_megacam_files execution started','MSC01_sort_megacam_files execution started.')
        
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'MSC01_sort_megacam_files')

        mcd.create_directory(base_path+'incomplete_files/',path_logfile,path_errorfile)
        
        os.chdir(base_path)
        for fits_filename in sorted(glob.glob('*.fits.fz')):
            date_obs = ''
            try:
                date_obs,keep_going = extract_date_megacam(fits_filename,keep_going)
            except:
                os.rename(fits_filename,base_path+'incomplete_files/'+fits_filename)
            if len(date_obs) == 10 and date_obs[4] == '-' and date_obs[7] == '-':
                dir_name = 'ut{:s}{:s}{:s}/'.format(date_obs[0:4],date_obs[5:7],date_obs[8:10])
                mcd.create_directory(base_path+dir_name,path_logfile,path_errorfile)
                mcd.output_log_entry(path_logfile,'Moving {:s} to {:s}...'.format(fits_filename,dir_name))
                os.rename(fits_filename,base_path+dir_name+fits_filename)

        mcd.send_status_email('MSC01_sort_megacam_files execution complete','MSC01_sort_megacam_files execution complete.')

    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')
    
    return None

if __name__ == '__main__':
    main()
