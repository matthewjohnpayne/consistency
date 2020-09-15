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
import copy

# ------------------------------------------------------------------------
# RAY PARALLELIZATION
# ------------------------------------------------------------------------
import ray
#try:
#    ray.init('auto')
#except:
#    ray.init()

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

    # Establish one-off connection to the database on mpcdb1
    # NB Despite being flat-file focused, we might need to
    # get submission IDs from the database ...
    host     = 'mpcdb1.cfa.harvard.edu'
    database = 'vmsops'
    cnx = psycopg2.connect(f"host={host} dbname={database} user=postgres")

    # Establish the list of objects / designations / numbers / ... to be queried
    desigs = [mc.unpacked_to_packed_desig(f'({x})') for x in range( n0,n1 )]

    # Check the consistency of each design
    #try:
    #    results =  ray.get(
    #                    [    establish_internal_consistency_of_flat_files_for_single_desig.remote( desig, DEBUG = True ) \
    #                    for desig in desigs ])
    #except:
    results =  [    establish_internal_consistency_of_flat_files_for_single_desig( desig, cnx, DEBUG = True ) \
                    for desig in desigs ]
                    
#@ray.remote
def establish_internal_consistency_of_flat_files_for_single_desig( desig, cnx, DEBUG = True ):
    '''
    '''
    # Get obs from flat files
    # - returns a list of obs
    print('Getting obs from ff for ', desig)
    obs_ff = ff.get_obs_from_ff(desig, DEBUG=DEBUG)

    # Look for duplicates
    print('Looking for any duplicates ')
    deduped_obs_list, duplicates = ff.find_duplicates(obs_ff)

    # Look for 2-line obs
    # Combine together where possible
    # Identify problems where exist
    print('Looking for orphan 2-line obs ...')
    obs_dict, orphans = ff.combine_two_line_obs(deduped_obs_list)

    # Are there other problems with the flat-file data that we can look for ?
    # (1) - Sometimes we do not have "Note 2" before 2020 in obs80: replace blank with default ?

    # Identify and fix any problems within the flat-file data
    if duplicates or orphans :
        
        # Fix simple duplicates ...
        if duplicates:
            print('Fixing duplicates ...')
            report = fix_primary_flat_file_data(desig, duplicates, [] , DELETING=True )
           
       
        # Attempt to fix sat/roving stuff here ...
        if orphans:
            print('Fixing orphans ...')
            correct_list, incorrect_list = [],[]
            for orphan in orphans:
            
                # Extract the obs80 bit ...
                obs80_bit       = orphan[15:56]
                print(f'orphan={orphan}, obs80_bit={obs80_bit}')
                
                # Find the original observations
                original_obs80_artifact  = ff.extract_original_observation(obs80_bit , ff.find_original_submission_artifact(obs80_bit, cnx) )
                print(f'original_obs80_artifact  : {original_obs80_artifact} ')
                print(f'incorrect_published_obs80  : {incorrect_published_obs80} ')
                
                # Construct a corrected obs80 bit
                corrected_obs80 = ff.construct_correct_obs80( \
                                        original_obs80_artifact ,
                                        incorrect_published_obs80 )
                print(f'corrected_obs80 : {corrected_obs80} ')

                correct_list.append(corrected_obs80)
                incorrect_list.append(incorrect_published_obs80)

            # Correct all of the orphans for a single flat file at once
            if incorrect_list != []:
                report = fix_primary_flat_file_data(desig, incorrect_list, correct_list )
                print(f' report from fix_primary_flat_file_data : {report} ')

            


# ------------------------------------------------------------------------
# MAKING DB CONSISTENT WITH FF
# - Older code that needs to be updated
# ------------------------------------------------------------------------

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
    
    
    # (i) Start by checking the internal self-consistency of the flat-files
    check_flat_file_internal_consistency(   int(sys.argv[1]),
                                            int(sys.argv[2]),
                                            DEBUG = True if len(sys.argv) > 3 and ( bool(sys.argv[3]) or sys.argv[3] == 'DEBUG') else False )
    
