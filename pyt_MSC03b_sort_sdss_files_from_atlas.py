import sys
import datetime
import glob, os, bz2, subprocess
import os.path
import macadamia_functions as mcd

def sort_file(file_tosort,log_file):
    data_dir             = int(file_tosort[8:14])
    data_subdir          = int(file_tosort[15:16])
    data_dir_path        = '/data1/data_sdss_processed/{:d}/'.format(data_dir)
    data_subdir_path     = '/data1/data_sdss_processed/{:d}/{:d}/'.format(data_dir,data_subdir)
    data_srctables_path  = '/data1/data_sdss_processed/{:d}/{:d}/source_tables/'.format(data_dir,data_subdir)
    data_photoutput_path = '/data1/data_sdss_processed/{:d}/{:d}/phot_output/'.format(data_dir,data_subdir)

    if not os.path.exists(data_dir_path):
        os.mkdir(data_dir_path)
    if not os.path.exists(data_subdir_path):
        os.mkdir(data_subdir_path)
    if not os.path.exists(data_srctables_path):
        os.mkdir(data_srctables_path)
    if not os.path.exists(data_photoutput_path):
        os.mkdir(data_photoutput_path)
        
    dir_destination = data_subdir_path

    if file_tosort[-8:] == '.fits.fz':
        dir_destination = data_subdir_path
    if file_tosort[-10:] == '.wcssolved':
        dir_destination = data_subdir_path
    if file_tosort[-11:] == '.photsolved':
        dir_destination = data_subdir_path
        
    if file_tosort[-26:] == '_multiap_photometry.txt.gz':
        dir_destination = data_srctables_path
    if file_tosort[-26:] == '_calibrated_sources.txt.gz':
        dir_destination = data_srctables_path
        
    if file_tosort[-11:] == '.moments.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-11:] == '.srclist.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-9:] == '.stars.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-10:] == '.trails.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-18:] == '_histogram1.pdf.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-18:] == '_histogram2.pdf.gz':
        dir_destination = data_photoutput_path
    if file_tosort[-22:] == '_refcat_sources.dat.gz':
        dir_destination = data_photoutput_path
    
    
    print('{:s} - Moving {:s} to {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),file_tosort,dir_destination))
    log_file.write('{:s} - Moving {:s} to {:s}...\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),file_tosort,dir_destination))
    mcd.move_file(file_tosort,dir_destination)
    return log_file


def main():

    # Define filenames and paths
    if len(sys.argv)!=2:
        print('Usage:\n python3 pyt_MSC03b_sort_sdss_files_from_atlas.py [staged_data_directory]\n')
    staged_data_directory = sys.argv[1]
    
    path_logfile = '/data1/log_files/log_{:s}_MSC03b_sort_sdss_files_from_atlas.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    os.chdir(staged_data_directory)
    with open(path_logfile,'w') as log_file:
        for dir_tosort in sorted(glob.glob('*/')):
            dir_tosort_path = staged_data_directory + dir_tosort
            os.chdir(dir_tosort_path)
            for file_tosort in sorted(glob.glob('frame-*')):
                log_file = sort_file(file_tosort,log_file)
            
    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
