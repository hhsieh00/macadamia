import re
import sys
sys.path.append('/Users/hhsieh/anaconda3/envs/astroconda/lib/python3.6/site-packages')
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages')
import glob, os
import datetime

def main():

    # Define filenames and paths
    if len(sys.argv)!=2:
        print('Usage:\n python3 pyt_MSC04_process_decam_ssois_results.py [base_path]\n')
    base_path = sys.argv[1]
    
    os.chdir(base_path)
    
    path_logfile = base_path + 'log_{:s}_MSC04_process_suprimecam_ssois_results.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    with open(path_logfile,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

    for ssois_filename in sorted(glob.glob('parsed_ssois_results_*_toproc.txt')):
        print('{:s} - Processing {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),ssois_filename))
        with open(path_logfile,'a') as log_file:
            log_file.write('{:s} - Processing {:s}...'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),ssois_filename))

        ssois_outputfilename = ssois_filename[:-11] + '_toingest.txt'

        with open(ssois_filename,'r') as ssois_file:
            with open(ssois_outputfilename,'w') as of:
                of.write('Object      Date       Time          Filter      Exptime  RA          Dec         Target                          Tel/Inst              Image Data Link\n')
            for _ in range(1): #skip first 1 header line
                next(ssois_file)
            for line in ssois_file:
                if line[37:43] == 'W-S-G+' or line[37:43] == 'W-S-R+' or line[37:43] == 'W-S-I+' or line[37:43] == 'W-S-Z+' or line[37:42] == 'W-J-B' or line[37:42] == 'W-J-V' or line[37:43] == 'W-C-RC' or line[37:43] == 'W-C-IC':
                    print('{:s} - {:s} ... {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),line[:31],line[375:394]))
                    with open(path_logfile,'a') as log_file:
                        log_file.write('{:s} - {:s} ... {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),line[:31],line[375:394]))
                    with open(ssois_outputfilename,'a') as of:
                        of.write(line)

    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
