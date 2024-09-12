import os, sys
import os.path
import glob, bz2, subprocess
import datetime
import math
import numpy as np
from astropy.io import fits
from astropy.io.fits import getheader
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import sqlite3
from sqlite3 import Error
import urllib
import smtplib, ssl


##### FUNCTION DEFINITIONS -- BASIC FILE FUNCTIONS #####

def create_directory(path,path_logfile,path_errorfile):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            output_error_log_entry(path_logfile,path_errorfile,'Creation of directory {:s} failed.'.format(path))
        else:
            output_log_entry(path_logfile,'Directory {:s} successfully created.'.format(path))
    return None

def rm_file(path_file):
    try:
        cmd_rm = '/bin/rm'
        cmd = [cmd_rm,path_file]
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: rm_file()')
        print(e)
    return None

def rename_file(old_filename,new_filename):
    try:
        cmd_mv = '/bin/mv'
        cmd = [cmd_mv,old_filename,new_filename]
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: rename_file')
        print(e)
    return None

def move_file(old_filename,new_filename):
    try:
        cmd_mv = '/bin/mv'
        cmd = [cmd_mv,old_filename,new_filename]
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: move_file')
        print(e)
    return None

def copy_file(old_filename,new_filepath):
    try:
        cmd_cp = '/bin/cp'
        cmd = [cmd_cp,old_filename,new_filepath]
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: copy_file')
        print(e)
    return None

def decompress_file_bzip2(path_file):
    try:
        bzip1,bzip2 = 'bzip2','-d'
        cmd = [bzip1,bzip2,path_file]
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: decompress_file_bzip2')
        print(e)
    return None

def compress_file_gzip(path_file):
    try:
        cmd_gzip = 'gzip'
        cmd = [cmd_gzip,path_file]
        if os.path.exists(path_file+'.gz'):
            os.rename(path_file+'.gz',path_file+'.backup.gz')
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: compress_file_gzip')
        print(e)
    return None

def decompress_file_gzip(path_file):
    try:
        cmd_gunzip = 'gunzip'
        cmd = [cmd_gunzip,path_file]
        if os.path.exists(path_file[:-3]) and path_file[-3:] == '.gz':
            os.rename(path_file[:-3],path_file[:-3]+'.backup')
        process = subprocess.call(cmd)
    except Exception as e:
        print('Function failed: decompress_file_gzip')
        print(e)
    return None

def compress_file_fpack(path_file):
    try:
        fpack_cmd = ''
        if os.path.isfile('/usr/local/bin/fpack'):
            fpack_cmd = '/usr/local/bin/fpack'
        elif os.path.isfile('/Users/hhsieh/Astro/tools/cfitsio/fpack'):
            fpack_cmd = '/Users/hhsieh/Astro/tools/cfitsio/fpack'
        elif os.path.isfile('/atlas/bin/fpack'):
            fpack_cmd = '/atlas/bin/fpack'
        if fpack_cmd != '':
            delete_cmd = '/bin/rm'
            cmd = [fpack_cmd,path_file]
            process = subprocess.call(cmd)
            cmd = [delete_cmd,path_file]
            process = subprocess.call(cmd)
        else:
            print('Function failed: compress_file_fpack() (command not found)')
    except Exception as e:
        print('Function failed: compress_file_fpack()')
        print(e)
    return None

def decompress_file_funpack(path_file):
    try:
        funpack_cmd = ''
        if os.path.isfile('/usr/local/bin/funpack'):
            funpack_cmd = '/usr/local/bin/funpack'
        elif os.path.isfile('/Users/hhsieh/Astro/tools/cfitsio/funpack'):
            funpack_cmd = '/Users/hhsieh/Astro/tools/cfitsio/funpack'
        elif os.path.isfile('/atlas/bin/funpack'):
            funpack_cmd = '/atlas/bin/funpack'
        if funpack_cmd != '':
            delete_cmd = '/bin/rm'
            cmd = [funpack_cmd,path_file]
            process = subprocess.call(cmd)
            cmd = [delete_cmd,path_file]
            process = subprocess.call(cmd)
        else:
            print('Function failed: compress_file_funpack() (command not found)')
    except Exception as e:
        print('Function failed: decompress_file_funpack()')
        print(e)
    return None

def sky2pix_fzfile(fits_file,ra,dec):
    xpix,ypix = -1,-1
    try:
        cmd_sky2pix = ''
        if os.path.isfile('/usr/local/bin/sky2pix'):
            cmd_sky2pix = '/usr/local/bin/sky2pix'
        elif os.path.isfile('/users/hhsieh/dropbox/hhsieh_db/astro_db/tools/wcseval/sky2pix'):
            cmd_sky2pix = '/users/hhsieh/dropbox/hhsieh_db/astro_db/tools/wcseval/sky2pix'
        elif os.path.isfile('/atlas/bin/sky2pix'):
            cmd_sky2pix = '/atlas/bin/sky2pix'
        if cmd_sky2pix != '':
            ra_str  = '{:.10f}'.format(ra)
            dec_str = '{:.10f}'.format(dec)
            cmd = [cmd_sky2pix,fits_file,ra_str,dec_str]
            process = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            out,err = process.communicate()
            output = out.split()
            xpix = float(output[0])
            ypix = float(output[1])
        else:
            print('Function failed: sky2pix_fzfile (command not found)')
    except Exception as e:
        print('Function failed: sky2pix_fzfile')
        print(e)
    return xpix,ypix


def compute_angular_dist_radec(ra1,dec1,ra2,dec2):
    # Compute angular distance in degrees from two sets of RA,Dec coordinates, also in degrees
    source_distance_deg = -1
    try:
        c1 = SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
        c2 = SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
        source_distance_deg = c1.separation(c2).degree
    except Exception as e:
        print('Function failed: compute_angular_dist()')
        print(e)
    return source_distance_deg


def get_pixel_coords_from_radec(fits_file,ra,dec):
    px,py = [-1,-1],[-1,-1]
    try:
        header = fits.getheader(fits_file)
        w = WCS(header)
        px,py = w.wcs_world2pix(ra,dec,1)
    except Exception as e:
        print('Function failed: get_pixel_coords_from_radec()')
        print(e)
    return px,py


def get_radec_from_pixel_coords(fits_file,xcoord,ycoord):
    ra,dec = -90,-90
    try:
        header = fits.getheader(fits_file)
        w = WCS(header)
        ra,dec = w.wcs_pix2world(xcoord,ycoord,1)
        keep_going = True
    except Exception as e:
        print('Function failed: get_radec_from_pixel_coords()')
        print(e)
        send_status_email('Function failed: mcd.get_radec_from_pixel_coords()','Function failed: mcd.get_radec_from_pixel_coords()')
        keep_going = False
    return ra,dec


##### FUNCTION DEFINITIONS -- DATABASE FUNCTIONS #####

def create_connection(db_file):
    # create a database connection to the SQLite database specified by the db_file
    # param db_file: database file
    # return: Connection object or None
    try:
        conn = sqlite3.connect(db_file)
        conn.text_factory = str
        return conn
    except Error as e:
        print('Function failed: create_connection()')
        print(e)
    return None

def add_apostrophe_escape(ast_name):
    # add escape sequence for apostrophes for use in SQL insertion command
    # param ast_name: original asteroid name with unescaped apostrophe(s)
    # return: ast_name_new: modified asteroid name with escaped apostrophe(s)
    num_apostrophes = 0
    name_length = len(ast_name)
    for idx in range(0,name_length):
        if ast_name[idx] == "'":
            num_apostrophes += 1
    ast_name_temp = [0 for idx in range(name_length+num_apostrophes)] # initialize new name string
    idx_new = 0
    for idx in range(0,name_length):
        if ast_name[idx] != "'":
            ast_name_temp[idx_new] = ast_name[idx]
        else:
            ast_name_temp[idx_new]   = "'"
            ast_name_temp[idx_new+1] = "'"
            idx_new += 1
        idx_new += 1
    ast_name_new = ''.join(ast_name_temp)
    return ast_name_new


##### FUNCTION DEFINITIONS -- EMAILS #####

def send_status_email(subject,message_text):
    port = 465  # For SSL
    smtp_server = "smtp.gmail.com"
    sender_email = "macadamia.asteroid.archive@gmail.com"  # Enter your address
    receiver_email = "hhsieh@gmail.com"  # Enter receiver address
    recipient = "Henry Hsieh <hhsieh@gmail.com>"
    #password = input("Type your password and press enter: ")
    password = "D0n'teatthecrabd1p!"
    context = ssl.create_default_context()
    message = '{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),message_text)
    msg = 'To: {:s}\nSubject: {:s}\n\n{:s}'.format(recipient,subject,message)
    #with smtplib.SMTP_SSL(smtp_server, port, context=context) as server_ssl:
    #    server_ssl.login(sender_email,password)
    #    server_ssl.sendmail(sender_email,receiver_email,msg)
    return None

def send_status_email_old(subject,message):
    try:
        port = 465  # For SSL
        smtp_server = "smtp.gmail.com"
        sender_email = "macadamia.asteroid.archive@gmail.com"  # Enter your address
        receiver_email = "hhsieh@gmail.com"  # Enter receiver address
        recipient = "Henry Hsieh <hhsieh@gmail.com>"
        #password = input("Type your password and press enter: ")
        password = "D0n'teatthecrabd1p!"
        context = ssl.create_default_context()
        msg = 'To: {:s}\nSubject: {:s}\n\n{:s}'.format(recipient,subject,message)
        with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
            server.login(sender_email,password)
            server.sendmail(sender_email,receiver_email,msg)
    except Exception:
        print('{:s} - Function failed: send_status_email()'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    return None


##### FUNCTION DEFINITIONS -- INITIALIZING FUNCTIONS #####

def validate_input_params(base_path,sqlite_file):
    keep_going = True
    if base_path[-1:] != '/': base_path = base_path + '/'
    if not os.path.isdir(base_path):
        print('Directory {:s} not found.'.format(base_path))
        keep_going = False
    if not os.path.exists(sqlite_file):
        print('Database {:s} not found.'.format(sqlite_file))
        keep_going = False
    return base_path,keep_going

def initialize_log_file(base_path,function_name):
    # Create and open log file
    dir_logfiles = base_path + 'log_files/'
    if not os.path.isdir(dir_logfiles): os.mkdir(dir_logfiles)
    path_logfile = dir_logfiles + 'log_{:s}_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),function_name)
    with open(path_logfile,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    return path_logfile

def initialize_log_error_file(base_path,function_name):
    # Create and open log file
    dir_logfiles = base_path + 'log_files/'
    if not os.path.isdir(dir_logfiles): os.mkdir(dir_logfiles)
    path_logfile   = dir_logfiles + 'log_{:s}_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),function_name)
    path_errorfile = dir_logfiles + 'errors_{:s}_{:s}.txt'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M%S'),function_name)
    with open(path_logfile,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    with open(path_errorfile,'w') as log_file:
        log_file.write('Archival Asteroid Photometry Error Log: {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
    print('Log file path:   {:s}'.format(path_logfile))
    print('Error file path: {:s}'.format(path_errorfile))
    return path_logfile,path_errorfile


def output_log_entry(path_logfile,log_entry):
    with open(path_logfile,'a') as log_file:
        print('{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
        log_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
    return None

def output_log_entry_nonstring(path_logfile,log_entry):
    with open(path_logfile,'a') as log_file:
        print('{:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(log_entry)
        log_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(log_entry,file=log_file)
    return None

def output_error_log_entry(path_logfile,path_errorfile,log_entry):
    with open(path_logfile,'a') as log_file, open(path_errorfile,'a') as error_file:
        print('{:s} - {:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
        log_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
        error_file.write('{:s} - {:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S'),log_entry))
    return None

def output_error_log_entry_nonstring(path_logfile,path_errorfile,log_entry):
    with open(path_logfile,'a') as log_file, open(path_errorfile,'a') as error_file:
        print('{:s}'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(log_entry)
        log_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(log_entry,file=log_file)
        error_file.write('{:s}\n'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(log_entry,file=error_file)
    return None

def remove_error_log_if_empty(path_logfile,path_errorfile):
    if os.path.exists(path_logfile) and os.path.exists(path_errorfile):
        if len(open(path_errorfile).readlines()) < 2:
            os.remove(path_errorfile)
            output_log_entry(path_logfile,'Removing empty error log file.')
    return None

def remove_output_file_if_empty(path_logfile,path_outputfile):
    if os.path.exists(path_logfile) and os.path.exists(path_outputfile):
        if len(open(path_outputfile).readlines()) < 2:
            os.remove(path_outputfile)
            output_log_entry(path_logfile,'Removing empty output file: {:s}'.format(path_outputfile))
    return None


##### FUNCTION DEFINITIONS -- WCS FUNCTIONS #####

def get_pixel_coords_from_radec(fits_file,ra,dec):
    px,py = [-1,-1],[-1,-1]
    try:
        header = fits.getheader(fits_file)
        w = WCS(header)
        px,py = w.wcs_world2pix(ra,dec,1)
    except Exception as e:
        print('Function failed: get_pixel_coords_from_radec()')
        print(e)
    return px,py


def get_radec_from_pixel_coords(fits_file,xcoord,ycoord):
    ra,dec = -90,-90
    try:
        header = fits.getheader(fits_file)
        w = WCS(header)
        ra,dec = w.wcs_pix2world(xcoord,ycoord,1)
        keep_going = True
    except Exception as e:
        print('Function failed: get_radec_from_pixel_coords()')
        print(e)
        keep_going = False
    return ra,dec


def compute_angular_dist(ra1,dec1,ra2,dec2):
    # Compute angular distance in degrees from two sets of RA,Dec coordinates, also in degrees
    source_distance_deg = -1
    try:
        c1 = SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
        c2 = SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
        source_distance_deg = c1.separation(c2).degree
    except Exception as e:
        print('Function failed: compute_angular_dist()')
        print(e)
    return source_distance_deg


def magavg(mag_array,magerr_array):
    avgmag,avgmagerr = -999,-999
    try:
        num_mags = len(mag_array)
        intensity_total = 0
        error_total = 0
        # initialize intensity array
        intensity = [[0 for idx1 in range(2)] for idx2 in range(num_mags)]
        for idx in range(0,num_mags):
            intensity[idx][0] = 10**(0.4*(0.-mag_array[idx]))
            intensity[idx][1] = ((magerr_array[idx]**2) * ((-0.921*(10**(-0.4*mag_array[idx])))**2.))**0.5
            intensity_total   = intensity_total + intensity[idx][0] / (intensity[idx][1]**2)
            error_total       = error_total + (intensity[idx][1]**2)**(-1)
        avgintensity = intensity_total / error_total
        avgintensityerr = error_total ** (-0.5)
        avgmag = -2.5 * math.log10(avgintensity)
        avgmagerr = 1.08574 * avgintensityerr/avgintensity
    except Exception as e:
        print('{:s} - Function failed: magavg()'.format(datetime.datetime.today().strftime('%Y-%m-%d %H:%M:%S')))
        print(e)
    return(avgmag,avgmagerr)


