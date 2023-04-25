import re
import sys
import datetime
import math
import urllib
import http.client
import glob, os, bz2, subprocess
import os.path
import requests
import smtplib,ssl
import macadamia_functions as mcd


##### FUNCTION DEFINITIONS -- DATABASE FUNCTIONS #####

def execute_parse_ssois_queries(search_urls_file,path_queryresults,parsed_results_directory,keep_going,path_logfile,path_errorfile):
    # Execute list of SSOIS queries
    if keep_going:
        try:
            mcd.output_log_entry(path_logfile,'Executing SSOIS queries listed in {:s}'.format(search_urls_file))
            parsed_ssois_results_filename = parsed_results_directory + 'parsed_ssois_results_' + search_urls_file[12:-4] + '_toingest.txt'
            if not os.path.isfile(parsed_ssois_results_filename):
                with open(parsed_ssois_results_filename,'w') as parsed_ssois_results_file:
                    parsed_ssois_results_file.write('Object      Date       Time          Filter      Exptime  RA          Dec         Target                          Tel/Inst              Image Data Link\n')
            with open(search_urls_file) as url_file:
                for _ in range(1): #skip first header line
                    next(url_file)
                for line in url_file:
                    obj_desig,epoch1,epoch2 = parse_ssois_url(line.strip(),path_logfile,path_errorfile)
                    ssois_search_result_filepath = create_ssois_result_path(path_queryresults,obj_desig,epoch1,epoch2)
                    search_executed_already = execute_ssois_query(line.strip(),ssois_search_result_filepath,path_logfile,path_errorfile)
                    if not search_executed_already:
                        parse_ssois_results_file(path_queryresults,ssois_search_result_filepath,parsed_ssois_results_filename,path_logfile,path_errorfile)
            mcd.compress_file_gzip(parsed_ssois_results_filename)
            mcd.output_log_entry(path_logfile,'SSOIS queries executed and results saved to {:s}'.format(parsed_ssois_results_filename+'.gz'))
        except Exception as e:
            mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_queries()'.format(search_urls_file))
            keep_going = False
            print(e)
    else:
        mcd.output_log_entry(path_logfile,'Function skipped: execute_ssois_queries()')
        keep_going = False
    return keep_going


def parse_ssois_url(ssois_url,path_logfile,path_errorfile):
    # Parse SSOIS URL
    obj_desig,epoch1,epoch2 = '','',''
    try:
        mcd.output_log_entry(path_logfile,'Parsing SSOIS URL {:s}...'.format(ssois_url))
        strlength = len(ssois_url)
        idx = 0
        while idx < strlength:
            if idx < (strlength-7):
                if ssois_url[idx:idx+7] == 'object=':
                    idx1 = idx+7
                    while ssois_url[idx1] != '&':
                        if ssois_url[idx1] != '+':
                            obj_desig = obj_desig + ssois_url[idx1]
                        idx1 += 1
            if idx < (strlength-7):
                if ssois_url[idx:idx+7] == 'epoch1=':
                    idx1 = idx+7
                    while ssois_url[idx1] != '&':
                        if ssois_url[idx1] != '+':
                            epoch1 = epoch1 + ssois_url[idx1]
                        idx1 += 1
            if idx < (strlength-7):
                if ssois_url[idx:idx+7] == 'epoch2=':
                    idx1 = idx+7
                    while ssois_url[idx1] != '&':
                        if ssois_url[idx1] != '+':
                            epoch2 = epoch2 + ssois_url[idx1]
                        idx1 += 1
            idx += 1
        mcd.output_log_entry(path_logfile,'Parsing SSOIS URL done.')
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed: parse_ssois_url() ({:s})'.format(ssois_url))
        print(e)
    return obj_desig,epoch1,epoch2


def create_ssois_result_path(path_queryresults,obj_desig,epoch1,epoch2):
    ssois_search_result_filepath = path_queryresults + 'ssois_results_{:07d}_{:s}_{:s}.txt'.format(int(obj_desig),epoch1,epoch2)
    return ssois_search_result_filepath


def extract_main_query_results_from_file(ssois_search_result_filepath,path_logfile,path_errorfile):
    mcd.output_log_entry(path_logfile,'Extracting main query results data from {:s}...'.format(ssois_search_result_filepath))
    objectname = ssois_search_result_filepath[-29:-22]
    with open(ssois_search_result_filepath) as df, open(ssois_search_result_filepath[:-4]+'_temp.txt','w') as of:
        found_start_main_content,found_end_main_content = False,False
        for line in df:
            if not found_start_main_content and len(line) > 43:
                if line[12:44] == '<!-- MAIN CONTENT begins here-->':
                    found_start_main_content = True
            if found_start_main_content and not found_end_main_content:
                of.write(line)
                if len(line) > 10:
                    if line[:10] == 'Query took':
                        found_end_main_content = True
    os.remove(ssois_search_result_filepath)
    os.rename(ssois_search_result_filepath[:-4]+'_temp.txt',ssois_search_result_filepath)
    mcd.output_log_entry(path_logfile,'Extracting main query results data from {:s} complete'.format(ssois_search_result_filepath))
    return None


def execute_ssois_query(url,result_filepath,path_logfile,path_errorfile):
    search_executed_already = False
    try:
        result_filepath_gz = result_filepath + '.gz'
        if not os.path.isfile(result_filepath_gz):
            mcd.output_log_entry(path_logfile,'Executing SSOIS query...{:s}'.format(url))
            query_successful = True
            max_retries = 50
            for _ in range(max_retries):
                try:
                    urllib.request.urlretrieve(url.strip(),result_filepath)
                except Exception:
                    mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS query failed. Retrying...'.format(url))
                else:
                    break
            else:
                mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS query failed ({:s}). Maximum retries reached'.format(url))
                query_successful = False
            if query_successful:
                extract_main_query_results_from_file(result_filepath,path_logfile,path_errorfile)
                mcd.output_log_entry(path_logfile,'Saving SSOIS query results to {:s}'.format(result_filepath))
                #mcd.compress_file_gzip(result_filepath)
        else:
            search_executed_already = True
    #except urllib.error.URLError as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    #except urllib.error.HTTPError as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    #except urllib.error.ContentTooShortError as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    #except http.client.IncompleteRead as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    #except http.client.HTTPException as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    #except ConnectionResetError as e:
    #    mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
    #    print(e)
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: execute_ssois_query()'.format(url))
        print(e)
    return search_executed_already


def parse_ssois_results_file(path_queryresults,ssois_search_result_filepath,parsed_ssois_results_filename,path_logfile,path_errorfile):
    try:
        mcd.output_log_entry(path_logfile,'Parsing {:s}...'.format(ssois_search_result_filepath))
        objectname = ssois_search_result_filepath[-29:-22]
        with open(ssois_search_result_filepath) as df:
            found_main_content,found_link,found_datetime,found_filter,found_exptime,found_ra,found_dec,found_target,found_telinst = False,False,False,False,False,False,False,False,False
            last_decam_file_string = ''
            no_images_found = False
            #rownum,newrow_rownum = 0,0
            for line in df:
                line_done = False
                #if len(line) > 43:
                #    if line[12:44] == '<!-- MAIN CONTENT begins here-->':
                if not found_start_main_content and len(line) > 50:
                    if line[0:30] == '<div class="table-responsive">':
                        found_main_content = True                        
                if found_main_content and not found_link and not line_done and not no_images_found:
                    if line[0:12] == '<td><a href=' or line[0:27] == '<td><a rel="external" href=':
                        image_link,char_idx = '',0
                        if line[0:12] == '<td><a href=':
                            while line[13+char_idx] != "'":
                                image_link = image_link + line[13+char_idx]
                                char_idx += 1
                            while line[13+char_idx+2:13+char_idx+6] != '</a>':
                                file_name = file_name + line[13+char_idx+2]
                                char_idx += 1
                        elif line[0:27] == '<td><a rel="external" href=':
                            while line[28+char_idx] != "'":
                                image_link = image_link + line[28+char_idx]
                                char_idx += 1
                            while line[28+char_idx+2:28+char_idx+6] != '</a>':
                                file_name = file_name + line[28+char_idx+2]
                                char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Image link found: {:s}'.format(image_link))
                        found_link = True
                        line_done = True
                if found_link and not found_datetime and not line_done and not no_images_found:
                    if len(line) > 8:
                        if line[0:4] == '<td>' and line[4:8].isnumeric():
                            date = line[4:14]
                            time = line[15:23]
                            mcd.output_log_entry(path_logfile,'   -- Date and time found: {:s} {:s}'.format(date,time))
                            found_datetime = True
                            line_done = True
                        elif len(line) > 35:
                            if line[31:35].isnumeric():
                                date = line[31:41]
                                time = line[42:50]
                                mcd.output_log_entry(path_logfile,'   -- Date and time found: {:s} {:s}'.format(date,time))
                                found_datetime = True
                                line_done = True                        
                if found_datetime and not found_filter and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        filter_name,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            filter_name = filter_name + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Filter found: {:s}'.format(filter_name))
                        found_filter = True
                        line_done = True
                if found_filter and not found_exptime and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        exptime,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            exptime = exptime + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Exposure time found: {:s}'.format(exptime))
                        found_exptime = True
                        line_done = True
                if found_exptime and not found_ra and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        rightasc,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            rightasc = rightasc + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Right ascension found: {:s}'.format(rightasc))
                        found_ra = True
                        line_done = True
                if found_ra and not found_dec and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        declination,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            declination = declination + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Declination found: {:s}'.format(declination))
                        found_dec = True
                        line_done = True
                if found_dec and not found_target and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        target,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            target = target + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Target found: {:s}'.format(target))
                        found_target = True
                        line_done = True
                if found_target and not found_telinst and not line_done and not no_images_found:
                    if line[0:4] == '<td>':
                        telinst,char_idx = '',0
                        while line[4+char_idx:9+char_idx] != '</td>':
                            telinst = telinst + line[4+char_idx]
                            char_idx += 1
                        mcd.output_log_entry(path_logfile,'   -- Telescope and instrument found: {:s}'.format(telinst))
                        found_telinst = True
                        line_done = True
                if found_telinst and line[0:5] == '</tr>' and not line_done and not no_images_found:
                    found_link,found_datetime,found_filter,found_exptime,found_ra,found_dec,found_target,found_telinst = False,False,False,False,False,False,False,False
                    if telinst == 'CTIO-4m/DECam':
                        if image_link[:33] == 'https://astroarchive.noirlab.edu/' and file_name[:4] == 'c4d_':
                            if file_name[18:21] == 'ooi':
                                decam_file_string = file_name[:24]
                                if decam_file_string != last_decam_file_string:
                                    with open(parsed_ssois_results_filename,'a') as parsed_ssois_results_file:
                                        #parsed_ssois_results_file.write('Object        Date       Time          JD              Tel/Inst              Exptime    Filter         Rdist     Ddist  PhsA  TrAnm  Vmag  Tmag  Nmag   Target                          PsAng  PsAMV  OrbPl  RA        Dec        RA          Dec          RA_rate  Dec_rate  RAsigma   DECsigma  EclLon  EclLat  GlxLon  GlxLat  Image_Data_Link\n')
                                        #parsed_ssois_results_file.write('{:>12s}  {:<10s} {:<12s}  {:14.6f}  {:<21s}  {:7.1f}  {:<10s}  {:8.3f}  {:8.3f}  {:4.1f}  {:5.1f}  {:4.1f}  {:4.1f}  {:4.1f}  {:<30s}  {:5.1f}  {:5.1f}  {:5.1f}  {:s}  {:s}  {:10.6f}  {:10.6f}  {:7.3f}  {:8.3f}  {:7.1f}  {:8.1f}  {:6.1f}  {:6.1f}  {:6.1f}  {:6.1f}  {:s} {:s}\n'.format(obj_desig,date,time,obs_datetime_mid.jd,telinst,float(exptime),filter_name,heliodist,geodist,phsang,trueanom,v_mag,t_mag,n_mag,target,sun_obj_pa,velocity_pa,orbpl_angle,ra_hms,dec_dms,float(rightasc),float(declination),ra_rate,dec_rate,ra_sigma,dec_sigma,ecllon,ecllat,glxlon,glxlat,file_name,image_link))
                                        parsed_ssois_results_file.write('{:>10s}  {:<10s} {:<12s}  {:<10s}  {:7.1f}  {:10.6f}  {:10.6f}  {:<30s}  {:<20s}  {:s} {:s}\n'.format(objectname,date,time,filter_name,float(exptime),float(rightasc),float(declination),target,telinst,file_name,image_link))
                                        mcd.output_log_entry(path_logfile,'{:>10s}  {:<10s} {:<12s}  {:<10s}  {:7.1f}  {:10.6f}  {:10.6f}  {:<30s}  {:<20s}  {:s} {:s}'.format(objectname,date,time,filter_name,float(exptime),float(rightasc),float(declination),target,telinst,file_name,image_link))
                                last_decam_file_string = decam_file_string
                    else:
                        with open(parsed_ssois_results_filename,'a') as parsed_ssois_results_file:
                            parsed_ssois_results_file.write('{:>10s}  {:<10s} {:<12s}  {:<10s}  {:7.1f}  {:10.6f}  {:10.6f}  {:<30s}  {:<20s}  {:s}\n'.format(objectname,date,time,filter_name,float(exptime),float(rightasc),float(declination),target,telinst,image_link))
                            mcd.output_log_entry(path_logfile,'{:>10s}  {:<10s} {:<12s}  {:<10s}  {:7.1f}  {:10.6f}  {:10.6f}  {:<30s}  {:<20s}  {:s}'.format(objectname,date,time,filter_name,float(exptime),float(rightasc),float(declination),target,telinst,image_link))
                    line_done = True
                    found_link,found_datetime,found_filter,found_exptime,found_ra,found_dec,found_target,found_telinst = False,False,False,False,False,False,False,False
                # Identify asteroids for which no matching images were found
                if found_main_content and len(line) > 28:
                    if line[0:29] == 'No matching images were found':
                        mcd.output_log_entry(path_logfile,'{:>10s}  -- No matching images found --'.format(objectname))
                        no_images_found = True
        mcd.compress_file_gzip(ssois_search_result_filepath)
    except Exception as e:
        mcd.output_error_log_entry(path_logfile,path_errorfile,'Function failed for {:s}: parse_ssois_results_file()'.format(ssois_search_result_filepath))
    return None


def main():

    # Define filenames and paths
    if len(sys.argv)!=4:
        print('Usage:\n python3 pyt_ingest_exposure_data_sdss.py [base_path] [sqlite_file] [thread_index]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path   = sys.argv[1]
    sqlite_file = sys.argv[2]
    thread_idx  = int(sys.argv[3])
    
    # Validate input parameters
    base_path,keep_going = mcd.validate_input_params(base_path,sqlite_file)

    if keep_going:
        mcd.send_status_email('SQ03_execute_parse_ssois_queries_{:02d} execution started'.format(thread_idx),'SQ03_execute_parse_ssois_queries_{:02d} execution started.'.format(thread_idx))
    
        # Create and initialize log file
        path_logfile,path_errorfile = mcd.initialize_log_error_file(base_path,'SQ03_execute_ssois_queries_{:02d}'.format(thread_idx))
    
        ssois_directory          = base_path + 'ssois_queries/'
        search_urls_directory    = ssois_directory + 'search_urls/'
        query_results_directory  = ssois_directory + 'query_results/'
        parsed_results_directory = ssois_directory + 'parsed_results/'
    
        if not os.path.isdir(ssois_directory):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS query directory {:s} not found'.format(ssois_directory))
            mcd.send_status_email('SQ03_execute_parse_ssois_queries_{:02d} execution failed'.format(thread_idx),'SQ03_execute_parse_ssois_queries_{:02d} execution failed - SSOIS query directory {:s} not found.'.format(thread_idx,ssois_directory))
            keep_going = False
        if not os.path.isdir(search_urls_directory):
            mcd.output_error_log_entry(path_logfile,path_errorfile,'SSOIS search URLs directory {:s} not found'.format(search_urls_directory))
            mcd.send_status_email('SQ03_execute_parse_ssois_queries_{:02d} execution failed'.format(thread_idx),'SQ03_execute_parse_ssois_queries_{:02d} execution failed - SSOIS search URLs directory {:s} not found.'.format(thread_idx,search_urls_directory))
            keep_going = False

    if keep_going:
        mcd.create_directory(query_results_directory,path_logfile,path_errorfile)
        mcd.create_directory(parsed_results_directory,path_logfile,path_errorfile)
        os.chdir(search_urls_directory)
        for search_urls_file_gz in sorted(glob.glob('search_urls_*_{:02d}.txt.gz'.format(thread_idx))):
            mcd.output_log_entry(path_logfile,'Executing SSOIS queries in {:s}...'.format(search_urls_file_gz))
            mcd.decompress_file_gzip(search_urls_file_gz)
            search_urls_file = search_urls_file_gz[:-3]
            path_queryresults = query_results_directory + 'query_results_' + search_urls_file[12:-4] + '/'
            mcd.create_directory(path_queryresults,path_logfile,path_errorfile)
            keep_going = execute_parse_ssois_queries(search_urls_file,path_queryresults,parsed_results_directory,keep_going,path_logfile,path_errorfile)
            if keep_going:
                search_urls_file_executed = search_urls_file[:-4]+'_executed.txt'
                os.rename(search_urls_file,search_urls_file_executed)
                mcd.compress_file_gzip(search_urls_file_executed)
            else:
                mcd.compress_file_gzip(search_urls_file)
                    
    mcd.send_status_email('SQ03_execute_parse_ssois_queries_{:02d} execution complete.'.format(thread_idx),'SQ03_execute_parse_ssois_queries_{:02d} execution complete.'.format(thread_idx))
    mcd.remove_error_log_if_empty(path_logfile,path_errorfile)
    mcd.output_log_entry(path_logfile,'Done.')

    return None


if __name__ == '__main__':
    main()

