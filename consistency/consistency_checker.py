#!/usr/bin/env python3
'''
Sketch-out of code to check (& fix/establish) observational consistency 

Items to discuss 
2020713 : 
(i) THINGS IN db BUT NOT IN ff :
 - Should we delete from db, or update the ff ? 
 - Eg (1) asterisk in line[13:14]
 - Eg (2) band in line[71:72] *** SEEMS TO BE PARTICULARLY A PROBLEM WITH T05 & T08 ***
    --- Specific Example can be seen for data from 02050 

20200714
Perhaps make a copy of the table (obs) structure 
Perhaps write new obs80 string
 - Would this be only the ones that changed, or everything ? 
Then populate derived fields
Then reconstruct obs80 from fields to demonstrate round-trip 
What would need to be frozen & when ? 
'''

# ------------------------------------------------------------------------
# THIRD PARTY IMPORTS
# ------------------------------------------------------------------------
import sys, os
import time
import psycopg2
import difflib
import numpy as np
from astropy.time import Time
from collections import defaultdict
import textwrap 
import healpy as hp
import glob
import re
import shutil
import random
import copy

# ------------------------------------------------------------------------
# RAY/DASK PARALLELIZATION
# ------------------------------------------------------------------------
#import ray
#ray.init('auto')
import dask
from distributed import Client
client = Client('tcp://131.142.192.121:8786')

# ------------------------------------------------------------------------
# MPC IMPORTS
# ------------------------------------------------------------------------
import mpc_psql as psql
import mpc_convert as mc
import mpc_new_processing_sub_directory as newsub
import status as mpc_status
sys.path.insert(0, '/share/apps/obs80')
import obs80 as o80

# ------------------------------------------------------------------------
# LOCAL IMPORTS
# ------------------------------------------------------------------------
import flat_files as ff

# ------------------------------------------------------------------------
# CROSS DESIGNATION DUPLICATES...
# ------------------------------------------------------------------------
   
def search_for_cross_designation_duplicates():
    '''
    There's a possibility that the same observation has
    been published against multiple object-designations
    Let's check for that
    
    NB I am not explicitly checking for duplicates WITHIN files here:
     - I am assuming I do that elsewhere ...
    '''
    print('search_for_cross_designation_duplicates')
    
    # We want to 'permanently' save some output files ...
    save_dir = '/sa/conchecks/data_products/'

    # ------------ NUMBERED FILES ------------------
    # Primary, published files
    files_ = glob.glob(f'/sa/mpn/N*dat', recursive=True)
    
    files_.extend(glob.glob(f'/sa/mpn/N*ADD', recursive=True))
    
    # In-progress ( between monthly pubs) files are in different location ...
    files_.extend( glob.glob(f'/sa/obs/*num', recursive=True) )
    
    # Save the num:file mapping, just in case ...
    file_dict = { n:f for n,f in enumerate(files_)}
    num       = { n:True for n,f in file_dict.items() } # Later on might want unnum files as well
    
    # ---------------- UN-numbered FILES -----------
    files_.extend(glob.glob(f'/sa/mpu/*dat', recursive=True))
    files_.extend(glob.glob(f'/sa/mpu/*ADD', recursive=True))
    
    file_dict = { n:f for n,f in enumerate(files_)}
    num       = { n:True if n in num else False for n,f in file_dict.items()}

    # ---------------- File-Mapping -----------
    filepath = os.path.join(save_dir,'file_num_mapping.txt')
    with open( filepath,'w') as fh:
        for n,f in file_dict.items():
            fh.write(f'{n},{f},{num[n]}\n')
    print('created...', filepath)
    
    
    
    '''
    # ---------------- Big data read ----------
    # Read the data into a single, massive dictionary
    # This is going to be challenging
    ALL = {}
    DUP = defaultdict(list)
    
    for n,f in file_dict.items():
        print(n,f)
        
        with open(f,'r') as fh:
            # local dict maps obs80-bit to integer representing file
            # NB: ignoring second-line obs, because those are the same for many many detections in the same exposure
            local     = {line[15:56]:n for line in fh if line[14] not in ['s','v']}
            
        # intersecn indicates duplicate obs80-bits
        intersecn = local.keys() & ALL.keys()
        
        # store duplicates with list of file-integers
        for k in intersecn:
            DUP[k].append(local[k])
            if isinstance(ALL[k], int):
                DUP[k].append(ALL[k])
            else:
                DUP[k].extend(ALL[k])

        # update the overall dictionary with local data
        ALL.update(local)
        
        # update the overall dictionary with the duplicates
        ALL.update(DUP)
        print('\t',len(ALL), len(DUP))

        # do a sanity print-out of the last input obs80bit
        lastkey = list(local.keys())[-1]
        print('\t'*2,' ...last key:value',lastkey, ALL[lastkey])
        
        del local
        del intersecn
        
        # Because I am impatient, I will print out the entire dict any time there is content ...
        if DUP:
            # save the duplicates to file
            filepath=os.path.join(save_dir,'duplicates.txt')
            with open( filepath,'w') as fh:
                for obs80bit, lst in DUP.items():
                    for i,n in enumerate(lst):
                        fh.write(f'{obs80bit},{i},{file_dict[n]},{num[n]}\n')
            print('\t'*3,'created/updated:', filepath)

        
    del ALL
    '''

    """
    # Get the list of pairs we need to check
    pairs = []
    for i in range(len(file_dict)):
        for j in range(i+1,len(file_dict)):
            pairs.append( (i,j) )
    print(f'len(file_dict)={len(file_dict)}')
    print(f'Need to check {len(pairs)} pairs to find duplicates...')
    random.shuffle(pairs)
    print('\t e.g.', pairs[:10])
    """
    
    # For each file, get the list of other files agsinst which it must be checked
    # NB "triangular", so don't repeat pair-wise comparisons.
    list_of_lists = [ (i,list(range(i+1,len(file_dict)))) for i in range(len(file_dict)) ]

    # ------------ FILE READ --------------------
    
    # Read all of the observations in a parallel style-ee
    for tup in list_of_lists:
        i   = tup[0]
        lst = tup[1]
        print('looping ... ',i, len(lst))
        
        # Read the i-file
        with open(file_dict[i],'r') as fh1:
            d1 = {line[15:56]:True for line in fh1 if line[14] not in ['s','v']}
            
        # Compare the i-file against all other files in the "lst" list
        # NB Processing in chunks
        list_of_dup_dicts_for_i =  []
        chunk = 200
        for k in range(0,len(lst),chunk):
            print(f'chunk...{k,len(lst),chunk}')
            chunk_lst = lst[k:k+chunk]
            list_of_dup_dicts_for_chunk = dask.compute( [ff.compare_file_against_provided_file_dict.remote(d1, file_dict[i2], i,i2) for i2 in chunk_lst])
            
            # Combine all presented dictionaries into a single dictionary for chunk
            list_of_dup_dicts_for_i.append( combine_list_dup_dicts(list_of_dup_dicts_for_chunk) )
        
        # Combine all presented dictionaries into a single dictionary
        list_of_dup_dicts.append( combine_list_dup_dicts(list_of_dup_dicts_for_i))
        
        print('\t-------',len(list_of_dup_dicts_for_i),len(list_of_dup_dicts[-1]))
        sys.exit()

    """
    list_of_dup_dicts = []
    chunk = 200
    for k in range(0,len(pairs),chunk):
    
        # For clarity, define which pairs we will look at in this chunk
        # NB : These are just pairs of integers
        chunk_pairs = pairs[k:k+chunk]
        print(f'chunking ... k={k}/{len(pairs)}, chunk={chunk}')
        
        #Get a list (of length chunk=200), where each entry is a dictionary of duplicates
        # - NB(1) dicts can be empty.
        # - NB(2) Parallelized over chunk=200 CPUs
        list_of_dup_dicts_for_chunk = dask.compute( [ff.check_two_files_for_dups.remote(    file_dict[p[0]],
                                                                                    file_dict[p[1]],
                                                                                    p[0],
                                                                                    p[1]) for p in chunk_pairs ] )
    
        # Combine all presented dictionaries into a single dictionary
        list_of_dup_dicts.append( combine_list_dup_dicts(list_of_dup_dicts_for_chunk))
        print('\t-------',len(list_of_dup_dicts[0]))
    """
        
    # Combine into a single dictionary
    DUPS = combine_list_dup_dicts(list_of_dup_dicts)
    
    # Save the duplicates to file
    if DUPS:
        print('\t'*3,'creating:', filepath)
        filepath=os.path.join(save_dir,'duplicates.txt')
        with open( filepath,'w') as fh:
            for obs80bit, lst in DUPS.items():
                for i,n in enumerate(lst):
                    fh.write(f'{obs80bit},{i},{file_dict[n]},{num[n]}\n')
        

    
def combine_list_dup_dicts(list_of_dup_dicts):
    # Combine all presented dictionaries into a single dictionary
    DUPS = defaultdict(list)
    for d in list_of_dup_dicts:
        for k,v in d.items():
            DUPS[k].extend(v) # <<-- v will always be a list
            
    # Ensure we don't have duplicates in DUPS
    return {k:list(set(v)) for k,v in DUPS.items()}

    
# ------------------------------------------------------------------------
# FLAT-FILE-ONLY CHECKS 
# ------------------------------------------------------------------------

def check_flat_file_internal_consistency( n0,n1, DEBUG = True ):
    '''
    High level function-call: used to check that flat-file observations are internally consistent

    inputs:
    -------
    n0,n1 : integers
     - first and last+1 "number" (permid) of numbered objects to be checked

    returns:
    -------
    nothing at present

    '''

    # Establish the list of objects / designations / numbers / ... to be queried
    desigs = [mc.unpacked_to_packed_desig(f'({x})') for x in range( n0,n1 )]
    
    # Establish a processing directory to work in ...
    proc_dir = newsub.generate_subdirectory( "obs_cons" )

    # Check the consistency of each design
    # (i) "ray" parallelization option
    '''
    if True:
        results =  dask.compute(
                        [    establish_internal_consistency_of_flat_files_for_single_desig.remote(  desig,
                                                                                                    proc_dir,
                                                                                                    cnx=None,
                                                                                                    DEBUG = True ) \
                        for desig in desigs ])
    #(ii) serial route
    else:
    
        # Establish one-off connection to the database on mpcdb1
        # NB Despite being flat-file focused, we might need to
        # get submission IDs from the database ...
        host     = 'mpcdb1.cfa.harvard.edu'
        database = 'vmsops'
        cnx = psycopg2.connect(f"host={host} dbname={database} user=postgres")

        #
        results =  [ establish_internal_consistency_of_flat_files_for_single_desig( desig,
                                                                                    proc_dir,
                                                                                    cnx=cnx,
                                                                                    DEBUG = True ) \
                        for desig in desigs ]
                    

    '''


# ------------------------------------------------------------------------
# MAKING DB CONSISTENT WITH FF
# - Older code that needs to be updated
# ------------------------------------------------------------------------
"""
def check_consistency(n0,n1, DEBUG = True):
    ''' High level function-call: used to check consistency of a list of numbered objects
    
    inputs:
    -------
    n0,n1 : integers
     - first and last+1 "number" (permid) of numbered objects to be checked

    returns:
    -------
    nothing at present

    '''
    
    # Establish one-off connection to the database on mpcdb1
    host = 'mpcdb1.cfa.harvard.edu'
    database = 'vmsops'
    cnx = psycopg2.connect(f"host={host} dbname={database} user=postgres")

    # Establish the list of objects / designations / numbers / ... to be queried 
    desigs = [mc.unpacked_to_packed_desig(f'({x})') for x in range( n0,n1 )]

    # Check the consistency of each design
    for desig in desigs:
        check_desig(cnx, desig, DEBUG=DEBUG) 
"""
"""
def check_desig(cnx, desig, DEBUG=False):
    ''' Declare a function to "process" / "fix" a *single* designation 

    inputs:
    -------
    cnx : psycopg2 connection object
    desig : unpacked number/permid of object being fixed 

    returns:
    --------
    nothing at present

    '''
    print('\n',desig)
    
    # Get obs from db
    # - returns a dict of obs, keyed on obs80_bit
    '''
    obs_db, probs_db = get_obs_from_db(cnx, desig, DEBUG=DEBUG)
    if probs_db:
        print_read_probs_dict(probs_db, 'DB')
    '''
    
    # Get obs from flat files 
    # - returns a dict of obs, keyed on obs80_bit
    obs_ff, probs_ff = ff.get_obs_from_ff(desig, DEBUG=DEBUG)
    
    
    # Identify and fix any problems within the flat-file data
    
    if probs_ff:
        print_read_probs_dict(probs_ff, 'FF')
        print()
        print(probs_ff['status'])
        print()
        # Fix simple duplicates ...
        lines_to_delete = []
        for status, incorrect_published_obs80 in zip(probs_ff['status'], probs_ff['line']):
            if status == -1 :
                lines_to_delete.append( incorrect_published_obs80 )
        if lines_to_delete != []:
            report = fix_primary_flat_file_data(desig, lines_to_delete, [] , DELETING=True )
            
        
        # For the sake of clarity while developing, will attempt to fix sat/roving stuff here ...
        # N.B. : "-2" is the code MJP is using to indicate a problem with 2-line data
        correct_list, incorrect_list = [],[]
        for status, incorrect_published_obs80 in zip(probs_ff['status'], probs_ff['line']):
            if status == -2 : 
                print('\nBefore finding original ... incorrect_published_obs80=', incorrect_published_obs80)
                obs80_bit       = incorrect_published_obs80[15:56]
                original_obs80_artifact  = extract_original_observation(obs80_bit , find_original_submission_artifact(obs80_bit, cnx) )
                print(f'original_obs80_artifact  : {original_obs80_artifact} ')
                print(f'incorrect_published_obs80  : {incorrect_published_obs80} ')
                corrected_obs80 = construct_correct_obs80( original_obs80_artifact , incorrect_published_obs80 )
                print(f'corrected_obs80 : {corrected_obs80} ')

                correct_list.append(corrected_obs80)
                incorrect_list.append(incorrect_published_obs80)

        # Correct all of the problems for a single flat file ...
        if incorrect_list != []:
            report = fix_primary_flat_file_data(desig, incorrect_list, correct_list )
            print(f' report from fix_primary_flat_file_data : {report} ')

    #sys.exit()
        
    # Compare observations
    # - returns nothing useful at present (empty dict) 
    '''
    result = compare_observations(cnx, obs_ff, obs_db, DEBUG=DEBUG)
    '''
"""

def print_read_probs_dict(probs_dict, header_str):
    status_dict = {-1:'DUP', -2:'S/R', -3:'MAL', -4:'???'}
    print(f'{header_str}:probs with read-observations..')
    for status, line in zip(probs_dict['status'], probs_dict['line']):
        print(f'\t{status_dict[status], line}')
        if status == -2:
            pass # sys.exit(' *** PROBLEM WITH SAT / ROVING DATA ***')


    
if __name__ == "__main__":
    #check_consistency(int(sys.argv[1]),
    #                  int(sys.argv[2]),
    #                  DEBUG = True if len(sys.argv) > 3 and ( bool(sys.argv[3]) or sys.argv[3] == 'DEBUG') else False)
    
    # Check the internal self-consistency of the flat-files
    #check_flat_file_internal_consistency(   int(sys.argv[1]),
    #                                        int(sys.argv[2]),
    #                                        DEBUG = True if len(sys.argv) > 3 and ( bool(sys.argv[3]) or sys.argv[3] == 'DEBUG') else False )
    
    # search_for_cross_designation_duplicates
    search_for_cross_designation_duplicates()
