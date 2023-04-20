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
    
    path_logfile = base_path + 'log_{:s}_MSC04_process_decam_ssois_results.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'))
    with open(path_logfile,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))

    for ssois_filename in sorted(glob.glob('parsed_ssois_results_*_toproc.txt')):
        print('Processing {:s}...'.format(ssois_filename))
        with open(path_logfile,'a') as log_file:
            log_file.write('Processing {:s}...'.format(ssois_filename))

        ssois_outputfilename = ssois_filename[:-11] + '_toingest.txt'

        prev_line = ''
        with open(ssois_filename,'r') as ssois_file:
            #of.write('Object      Date       Time          Filter      Exptime  RA          Dec         Target                          Tel/Inst              Image Data Link\\n')
            #for _ in range(1): #skip first 2 header lines
            #    next(ssois_file)
            for line in ssois_file:
                if prev_line == '':
                    with open(ssois_outputfilename,'w') as of:
                        of.write(line)
                    print(line[:-1])
                else:
                    if len(prev_line) > 210:
                        if prev_line[37:39] != 'VR' and prev_line[37:39] != 'Y ' and prev_line[196:199] == 'ooi' and line[178:201] != prev_line[178:201]:
                            print('{:s} - {:s} ... {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),prev_line[:31],prev_line[178:-1]))
                            with open(path_logfile,'a') as log_file:
                                log_file.write('{:s} - {:s} ... {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),prev_line[:31],prev_line[178:-1]))
                            with open(ssois_outputfilename,'a') as of:
                                of.write(prev_line)
                            prev_line_written = True
                        else:
                            prev_line_written = False
                prev_line = line
            if prev_line_written and prev_line[37:39] != 'VR' and prev_line[37:39] != 'Y ' and prev_line[196:199] == 'ooi':
                print('{:s} - {:s} ... {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),prev_line[:31],prev_line[178:-1]))
                with open(path_logfile,'a') as log_file:
                    log_file.write('{:s} - {:s} ... {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),prev_line[:31],prev_line[178:-1]))
                with open(ssois_outputfilename,'a') as of:
                    of.write(prev_line)

    with open(path_logfile,'a') as log_file:
        print('{:s} - Done.'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        log_file.write('{:s} - Done.\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
            
    return None


if __name__ == '__main__':
    main()
