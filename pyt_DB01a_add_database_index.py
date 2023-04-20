import re
import sys
import time
import datetime
import glob, os, bz2, subprocess
import csv
import os.path
import sqlite3
from sqlite3 import Error
import macadamia_functions as mcd


def main():
    
    # Define filenames and paths
    if len(sys.argv)!=3:
        print('Usage:\n python3 pyt_add_telinsts.py [base_path] [sqlite_file]\n')
        print(" (Trailing '/' needed in path specification)\n")
        exit()
    base_path = sys.argv[1]
    sqlite_file = base_path + sys.argv[2]
    
    #conn = mcd.create_connection(sqlite_file)  # Open connection to database file
    #cursor = conn.cursor()
    #add_index_query = "CREATE UNIQUE INDEX idx_proc_data_file ON exposures(proc_data_file)"
    #print(add_index_query)
    #cursor.execute(add_index_query)
    #conn.commit() # Commit changes
    #conn.close()  # Close connection to database file
    
    conn = mcd.create_connection(sqlite_file)  # Open connection to database file
    cursor = conn.cursor()
    add_index_query = "CREATE UNIQUE INDEX idx_mosaic_id_base_filename ON exposures(mosaic_element_id,base_filename)"
    cursor.execute(add_index_query)
    conn.commit() # Commit changes
    conn.close()  # Close connection to database file

    return None


if __name__ == '__main__':
    main()

