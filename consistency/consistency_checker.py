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
import difflib
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
# MPC IMPORTS
# ------------------------------------------------------------------------
import mpc_psql as psql
import mpc_convert as mcimport mpc_new_processing_sub_directory as newsub
import status as mpc_status
sys.path.insert(0, '/share/apps/obs80')
import obs80 as o80

# ------------------------------------------------------------------------
# LOCAL IMPORTS
# ------------------------------------------------------------------------
import flat_files as ff

# ------------------------------------------------------------------------
# TOP-LEVEL FUNCTIONS TO HANDLE THE OVERALL RUNNING OF CONSISTENCY CHECKS
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
    host = 'mpcdb1.cfa.harvard.edu'
    database = 'vmsops'
    cnx = psycopg2.connect(f"host={host} dbname={database} user=postgres")

    # Establish the list of objects / designations / numbers / ... to be queried
    desigs = [mc.unpacked_to_packed_desig(f'({x})') for x in range( n0,n1 )]

    # Check the consistency of each design
    for desig in desigs:
    
        # Get obs from flat files
        # - returns a dict of obs, keyed on obs80_bit
        obs_ff, probs_ff = ff.get_obs_from_ff(desig, DEBUG=DEBUG)


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


# ------------------------------------------------------------------------
# FUNCTIONS RELATED TO READING FROM FLAT-FILES
# ------------------------------------------------------------------------
def get_obs_from_ff(desig, DEBUG=False):
    '''
    # Get obs from flat files 
    # - Do we need to search using findn & findu on multiple "identifications" ?  
    # - Would be good to double-check nothing extra is returned 
    # - What files do / don't get searched?
    '''
    if DEBUG:
        print('ff-query is basic: needs to be expanded to use findu & other (linked) designations')
        
    # Query stolen from Dave (see 'getobs_num.py)
    obs = []
    outfile = desig + "_ff.obs"
    cmd = f"""$findn -o={outfile} -O={desig}"""
    os.system(cmd)
    if os.path.isfile(outfile):       
        with open(outfile) as f:
            lines = f.readlines()
            obs_list = [_.strip() for _ in lines]
    
        
    # Find any duplicate lines
    obs, probs,     = {}, defaultdict(list), []
    for n,line in enumerate(obs_list):
        obs80_bit = line[15:56]
        
        # put duplicates into problem-dict for now
        if obs80_bit in obs:
            probs['status'].append( -1 )
            probs['line'].append( line )
        else:
            obs[obs80_bit]=line
    
    # Now Find and combine 2-line obs ...
    # ... make dict key-ed on obs80 bit
    # NB : The check is being made under the
    # assumption that any dups have been removed ...
    obs = {} 
    for n,line in enumerate(deduped_obs_list):
        # Satellites & Roving
        if line[14] in ['S','V'] and obs_list[n+1][14] == line[14].lower()  :
            obs[obs80_bit] = line+obs_list[n+1]
        elif line[14] in ['s','v'] and obs_list[n-1][14] == line[14].upper()  :
            pass # because of previous logic
        elif line[14] in ['S','V'] and obs_list[n+1][14] != line[14].lower()  or \
             line[14] in ['s','v'] and obs_list[n-1][14] != line[14].upper()  :
            probs['status'].append( -2 )
            probs['line'].append( line )
            
        # Standard lines 
        else:
            obs[obs80_bit] = line
            
    return obs, probs

def get_original_submission_id(obs80_bit, cnx):
    ''' From obs80_bit string(s), get original submission_id from database (where possible) '''
    
    # Get submission IDs from db (looks like '2020-06-23T00:43:11.000_0000DdpR') 
    query = f"SELECT submission_id FROM obs WHERE obs80_bit='{obs80_bit}'"
    rows = psql.psql_execute_query(query, cnx)

    # We would ideally like there to be only a single obs80_bit match in the db ...
    if len(rows) == 1 :
        submission_id = rows[0][0]
        return submission_id
    elif len(rows) > 1 :
        sys.exit('Multiple submission ids returned for obs80_bit = %s : %r' % (obs80_bit, rows ) )        
    else:
        sys.exit('Could not find submission id for obs80_bit = %s' % obs80_bit)

        
def find_original_submission_artifact(obs80_bit, cnx):
    ''' From obs80_bit string(s), get original submission data (where possible) 

    NB : I wanted to use formatobs.txt, but sometimes doesn't exist, so using "artifact" instead
    # - 3317 does not have a "formatobs" file ... gaps in artifact file ... how did it get processed ?

    inputs:
    -------

    returns:
    --------

    '''
    # Get submission IDs from db (looks like '2020-06-23T00:43:11.000_0000DdpR')
    submission_id = get_original_submission_id(obs80_bit, cnx)
    splt=re.split('-|T', submission_id)
    yr,mn,dy = splt[0],splt[1],splt[2]
    files_ = glob.glob(f'/sa/archive/*/{yr}/{mn}/{dy}/{submission_id}/artifact', recursive=True)
    assert files_, f'no archive file found at /sa/archive/*/{yr}/{mn}/{dy}/{submission_id}/artifact'
        
    # Extract original data
    with open(files_[0],'r' ) as fh:
        data = fh.readlines()
    return data

def extract_original_observation( obs80_bit, artifact_data ):
    ''' Find the obs80_bit in the artifact_data, return all associated lines (i.e. 2 lines if satellite /roving) '''
    matches_ = [i for i,_ in enumerate(artifact_data) if obs80_bit in _ ]
    if len(matches_) == 1 :
        i = matches_[0]
        if artifact_data[i][14] in ['S','V']:
            obs80 = artifact_data[i:i+2]
        elif artifact_data[i][14] in ['s','v']:
            obs80 = artifact_data[i-1:i+1]
        else:
            obs80 = artifact_data[i:i+1]
    elif len(matches_) == 0 :
        sys.exit('Zero returns from artifact-match:%r' % obs80_bit)
    else:
        sys.exit('Multiple returns from artifact-match:%r' % obs80_bit)
    return obs80

def construct_correct_obs80( original_obs80 , incorrect_published_obs80 ):
    '''  Create a corrected version of published obs80 string(s) 

    The assumption is that the incorrect_published_obs80 string 
    has been used to find the original submissiont-artifact 
    (via *find_original_submission_artifact()* ), and the original_obs80
    extracted (using *extract_original_observation()* ).

    We now want to create a corrected version that 
    (a) has the correct 2-line sat/roving data
    (b) has the correct publishing info, etc.

    *** Can this use / be combined with the *categorize_differences()* function [see below] ??? ***

    inputs:
    -------
    original_obs80 :
    incorrect_published_obs80 : 

    returns:
    --------
    correct_obs80 : 


    '''
    # If length of original ==2, then satellite ...
    if len(original_obs80):

        # stitch components together
        original_obs80 = ''.join([_.strip('\n') for _ in original_obs80])
        print('joined original_obs80        =', original_obs80)
        # If missing a line, get second line from original
        if len(incorrect_published_obs80) == 80:
            correct_obs80 = incorrect_published_obs80 + original_obs80[80:]

        # Ensure that the two lines of the just-constructed "correct_obs80"...
        # ... have consistent pub-refs & obscodes
        correct_obs80 = correct_obs80[:152] + correct_obs80[72:80]

    # If initial version of corrrect_obs80 doesn't exist, then crap out as can't fix ... 
    try:
        print('extended___________________ ', correct_obs80)
    except:
        sys.exit('cannot fix ...', original_obs80 , incorrect_published_obs80)

    # Using the cat_diff routine that was original developed to compare flat-0files / db
    # NB Needs equal lengths ...
    R = categorize_differences( correct_obs80, original_obs80 , DEBUG = False)

    # If acceptable status returned, then use returned action as replace
    if R.status == 1:
        correct_obs80 = R.action
        print('corrected___________________ ', correct_obs80)
    else:
        R.pprint()
        print()
        sys.exit(f'R.status = {R.status} , R.info={R.info}')

    return correct_obs80 

def fix_primary_flat_file_data( desig, incorrect_list, correct_list, DELETING=False):
    '''

    *** Need to be really careful about this ***
    *** Need to do something like ... ***
    *** (i) Find the relevant primary data file [can be in /sa/mpn or in tot*] 
    *** (ii) Freeze the system (lock status)
    *** (iii) Copy primary data file to temp location (esp. while developing) 
    *** (iv) Find the location of the incorrect data in the primary data file 
    *** (v) Replace the incorrect data with the correct data (non-trivial : needs to have pubn-record, etc) 
    *** (vi) Do some sense checks of the difference between the initial and fixed versions 
    *** (vii)  write the data to the temp file 
    *** (viii) Replace the primary data with the fixed copy
    *** (ix) Unlock the system 

    inputs:
    -------

    returns:
    --------

    '''
    # We want to 'permanently' save some output files ...
    save_dir = '/sa/conchecks/data_products/'
    
    
    #*** (i) Find the relevant primary data file [can be in /sa/mpn or in tot*]
    src_files = []
    for incorrect_published_obs80 in incorrect_list:
        src_files.extend( find_primary_data_file( desig, incorrect_published_obs80 ) )
        src_files = list(set(src_files))
    print('src_files = ', src_files)
    assert len(src_files), f'No src_file could be found that contains the incorrect data ... incorrect_list={incorrect_list}'
    
    #*** (ii) Freeze the system (lock status)
    # ~~~~~~~~~~~~~ IF WE CRAP OUT AT ANY POINT BELOW WE NEED TO RELEASE THE LOCK ~~~~~~~~~~~~~~~~~~
    print('Setting mpc_temp_status')
    mpc_status.set_status("mpc_temp_status","MJP_FIXING_PRIMARY_FLAT_FILES")
    
    try:
        
        #*** (iii) Copy primary data file to temp location (esp. while developing)
        # I am allowing for the possibility that there are multiple files to be fixed ...
        dst_dir = newsub.generate_subdirectory( "obs_cons" ) 
        for src_file in src_files:
            dst_file= os.path.join(  dst_dir , os.path.basename(src_file) ); print('dst_file = ',dst_file)
            shutil.copyfile(src_file, dst_file )

            # Read the primary data file
            with open(dst_file,'r') as fh : data = fh.readlines()

            # Files to write to so that MR can update mysql
            bad_filepath  = os.path.join(save_dir, desig + '_bad.dat')
            good_filepath = os.path.join(save_dir, desig + '_good.dat')
            print(f' bad_filepath= {bad_filepath} , good_filepath= {good_filepath} ')
            with open(bad_filepath, 'w') as bad_fh:
                with open(good_filepath, 'w') as good_fh:

                    # If we are deleting duplicates ...
                    if DELETING :
                        
                        seen = {} 
                        fixed_data = []
                        incorrect_dict = {_:True for _ in incorrect_list}
                        for line in data:
                            # If the lines are to be deleted, record to tell MR so that the mysql database can be updated
                            if line.strip('\n') in incorrect_dict and line not in seen:
                                bad_fh.write(line)
                            # If we are keeping the line ...
                            else:
                                fixed_data.append(line)
                            # Record that we have seen the line so that we can stop ourselves deleting it twice!
                            seen[line]=True
                            
                            
                    # If not deleting, but doing replacement ... 
                    else: 
                        for incorrect_published_obs80, corrected_obs80 in zip(incorrect_list, correct_list) :
                            fixed_data = []

                            # Check the inputs 
                            assert corrected_obs80 not in ['', ' ', [], [''], ['','']], \
                                'corrected_obs80 = {corrected_obs80} : not sure that this routine can cope with such input ...'
                            assert isinstance(incorrect_published_obs80, str), f'incorrect_published_obs80 is of type {type(incorrect_published_obs80)}, rather than a string'
    
                            #*** (iv) Find the location of the incorrect data in the primary data file
                            line_num = [i for i,line in enumerate(data) if incorrect_published_obs80.strip() in line]
                            assert len(line_num) < 3, 'len(line_num)={len(line_num)} which is >=3 which seems like a suspiciously large number so I am terminating...'
        
                            #*** (v) Replace the incorrect data with the correct data (the correct data has been created earlier)
                            #        At the same time we also output the incorrect & correct data to some files to be used to update the MYSQL database
                            for n,line in enumerate(data):
                                if n not in line_num :
                                    # We keep the normal stuff as-is
                                    fixed_data.append(line)
                                else:
                                    # For removal from mysql
                                    bad_fh.write(line)

                                    if isinstance(corrected_obs80, str):
                                        l = corrected_obs80 if corrected_obs80[-1]=='\n' else corrected_obs80+'\n'
                                        # Corrected data for flat files 
                                        fixed_data.append(l)
                                        # Corrected data for mysql
                                        good_fh.write(l)
                                    elif isinstance(corrected_obs80, list):
                                        for _ in corrected_obs80:
                                            l = _ if _[-1]=='\n' else _+'\n'
                                            # Corrected data for flat files 
                                            fixed_data.append(l)
                                            # Corrected data for mysql
                                            good_fh.write(l)
                                    else:
                                        sys.exit(f'corrected_obs80 is of type{type(corrected_obs80)}: do not know how to process')

                            #*** (vi) Do some sense checks of the difference between the initial and fixed versions
                            assert len(fixed_data) - len(data) == len(line_num), 'Lengths do not make sense: {len(fixed_data),len(data),len(line_num)} '

                            # copy fixed data into data ready for next loop around ...
                            data = copy.deepcopy(fixed_data)
                        
 
            #*** (vii)  write the data to the temp file 
            replace_file = dst_file + 'replace'
            assert not os.path.isfile(replace_file), 'replacement file {replace_file} already exists which is bad'
            with open(replace_file,'w') as fh :
                for line in fixed_data:
                    l = line if line[-1]=='\n' else line+'\n'
                    fh.write(l)
            assert os.path.isfile(replace_file), 'replacement file {replace_file} does NOT exist which is bad'
        
            #*** (viii) Replace the primary data with the fixed copy
            print(f'replacing file={src_file} with file {replace_file} ')
            #shutil.copyfile(replace_file, src_file)
            
            #*** (ix) Recreate the index files if necessary
            #         NEED TO BE CAREFUL ABOUT THIS ...
            #         (a) Mike/Dave indicated this is only necessary if the file being altered is one of the permanent,
            #             master files, rather than one of the temp *tot* files
            #         (b) However, my inspection of /share/apps/mpec/publish_dou_mpec.sh, /share/apps/com/indexed/update.sh
            #             (and sub-scripts) suggests that there *ARE* some form of index files for the temp/pending/within-month files
            #         (c) TO gain some understanding, the monthly-prep rebuilds are done here : /sa/com/indexed/update.sh [SAME AS ABOVE]
            #         (d) Given that ... calls /share/apps/com/indexed/buildnumupd.sh

            #*** (xi) Remove / Tidy-up the temp files & temp dir
            #shutil.rmtree(dst_dir)
            #assert not os.path.isdir(dst_dir), f'dst_dir={dst_dir} still exists when it should NOT'


        #*** (ix) Unlock the system
        print('Unsetting the mpc_temp_status')
        mpc_status.set_status("mpc_temp_status","")

    
    except Exception as e: 
        print('\n'*2)
        print('EXCEPTION IN fix_primary_flat_file_data')
        print('\n'*2)
        print(e)
        print('\n'*2)
        print('Unsetting the mpc_temp_status as part of the EXCEPTION handling')
        mpc_status.set_status("mpc_temp_status","")

    return True

    
    
def find_primary_data_file( desig, incorrect_published_obs80 ):
    ''' #*** (i) Find the relevant primary data file [can be in /sa/mpn or in tot*] '''
    potential_source_files = []
    
    # Numbered master files look like N0123456.dat
    files_ = glob.glob(f'/sa/mpn/N*dat', recursive=True)
    files_.sort()
    num_str  = [os.path.basename(_)[1:-4] for _ in files_]
    try:
        desig = int(desig)
    except:
        desig = int(mc.convert_packed_perm_desig(desig))
    print( 'finding primary data file for ', desig )
    potential_source_files.append( [files_[i] for i,n in enumerate(num_str) if int(n)<=desig][-1] )
    
    # In-progress ( between monthly pubs) files are in different location ...
    potential_source_files.extend( glob.glob(f'/sa/obs/*num', recursive=True) )
    
    # Check which file(s) incorrect_published_obs80 is in
    src_files = []
    for f in potential_source_files:
        with open(f, "r") as fh:
            for line in fh:
                if incorrect_published_obs80.strip() in line:
                    src_files.append(f)
                    
    # NB I am deliberately *NOT* returning the line number here ...
    # ... because there is a TINY chance the file will be altered ...
    # ... between now and when I lock the file 
    return list(set(src_files))
    
# ------------------------------------------------------------------------
# FUNCTIONS RELATED TO READING FROM DATABASE
# ------------------------------------------------------------------------

def get_obs_from_db(cnx, desig, DEBUG=False):
    '''
    # Get obs from db	
    # - Do we need to search on multiple "identifications", perhaps in multiple fields ? 
    # ... perhaps do this as an option / only if standard search is "not good" 
    # - Would be good to be consistent 

    inputs:
    -------

    returns:
    --------

    '''
    if DEBUG:
        print('db-query is basic: needs to be expanded to look at other fields & other (linked) designations')

    # Query stolen from Dave (see 'getobs_num.py)
    qry = psql.get_obs_by_status(desig, 'Pp', cnx)

    # Make dict key-ed on obs80 bit
    obs, probs = {}, defaultdict(list)
    for line in qry :

        obs80_bit = line[15:56]

        try :
            Opt(line)
            PARSED = True
        except:
            PARSED = False
            
        # put any duplicates into prob-dict for now
        if obs80_bit in obs:
            probs['status'].append( -1 )
            probs['line'].append( line )            
        
        # Look for malformed satellite / roving
        # It would be bad if we had the start of a satellite line, but the end was missing/incomplete
        # It would also be bad if we had the end of a satellite line as the start ...
        elif line[14] in ['S','V'] and (len(line) < 160 or line[94] not in ['s','v'] )or \
             line[14] in ['s','v']:
            probs['status'].append(-2)
            probs['line'].append(line)

        elif not PARSED:
            probs['status'].append(-3)
            probs['line'].append(line)
            

        # Standard lines
        else:
            obs[obs80_bit] = line        

        
    return obs, probs

# ------------------------------------------------------------------------
# FUNCTIONS RELATED TO COMPARING OBS ACROSS DATABASE & FLAT-FILES
# ------------------------------------------------------------------------
# Experimenting with class to hold repairs to obs80
# Wanted to use dataclasses, but only py3.6 on marsden 
class repair(object):
    '''

    '''
    __slots__=['status', 'action', 'info'] 
    def __init__(self, ):
        self.status = []
        self.action = []
        self.info   = []
        
    def app(self, status, action, info):
        self.status.append(status)
        self.action.append(action)
        self.info.append(info)
        
    def summarize(self,):
        ''' May have multiple components. Boil-down to single status ... '''
        
        # If there are any untreated problems, set STATUS == -ve 
        # If anything was fixed, set +ve, otherwise = 0 if exactly OK
        # *** The above means we do NOTHING if there are ANY NEGATIVE ***
        self.status = np.min(self.status) if np.any( np.array(self.status) < 0 ) else np.max(self.status)

        # If there is anything we can fix, set the action ...
        if self.status > 0 :
            self.action = [_ for _ in self.action if _ is not None]

        # If there is stuff we can't fix, set the info
        if self.status < 0 :
            self.info = [_ for _ in self.info if _ is not None]

    def pprint(self,):
        self.summarize()
        print(self.status)
        if self.status > 0:
            try:
                for _ in self.action : print(textwrap.fill(_, width=85))
            except:
                for _ in self.action : print(_)
        if self.status < 0 :
            try:
                for _ in self.info : print(textwrap.fill(_, width=85))
            except:
                for _ in self.info : print(_)

                
def compare_observations(cnx, obs_ff, obs_db, DEBUG=False):
    ''' compare a set of observations 

    inputs:
    ------
    cnx - connection object to postgres db

    obs_ff - dictionary of observations from flat-files 
     - keyed on obs80_bit, values = obs80 (but with 2-lines joined into single line where appropriate)

    '''
    
    # (a) Compare the "obs80_bit" sections of each observation to establish matching observations : 
    #     Establish the number of objects that are : 
    #     db_only [only in database], 
    #     ff_only [only in flat-files], 
    #     db_and_ff [present in both flat-files and database]
    #     (N.B. ...  N_db = N__db_only + N_db_and_ff & N_ff = N_ff_only + N_db_and_ff)
    ff_keys   = set(obs_ff.keys())
    db_keys   = set(obs_db.keys())
    db_and_ff = ff_keys.intersection(db_keys)
    ff_only   = ff_keys - db_and_ff
    db_only   = db_keys - db_and_ff
    if DEBUG:
        print(f'compare_observations : ...')
        print(f' N_ff_keys  ={len(ff_keys)}')
        print(f' N_db_keys  ={len(db_keys)}')
        print(f' N_db_and_ff={len(db_and_ff)}')
        print(f' N_ff_only  ={len(ff_only)}')
        print(f' N_db_only  ={len(db_only)}')
        
    # (b) ***ff_only*** observations
    #     Search the db based on obs80_bit to see whether these can be recovered
    #     Are they published?
    #     Are they listed under a different designation?
    #     What action/query would need to be taken to "fix" the database ?
    #     - Likely to not just be obs80 string, but also other fields
    for obs80_bit in ff_only:
        R = process_ff_only(cnx, obs80_bit, obs_ff, DEBUG=DEBUG)
        if DEBUG and R.status != 1:
            R.pprint()
        
    
    
    # (c) ***db_and_ff*** observations [in both db & ff], 
    #     Which observation-strings are identical ? 
    #     For non-identical, which sections of the strings differ?
    #        (e.g. the name/provid strings, the reference strings, asterisks),
    #        hence what action/query would need to be taken to "fix" the database ?
    #     Are there other db fields (besides obs80 string) that also need fixing?
    for obs80_bit in db_and_ff:
        if obs_ff[obs80_bit] != obs_db[obs80_bit]:
            R = categorize_differences(obs_ff[obs80_bit] , obs_db[obs80_bit], DEBUG = DEBUG)
            if DEBUG and R.status != 1 :
                R.pprint()
                    
    # (d) For the N_db_only observations: what the hell do we do about these? 
    #     (1) Start by doing nothing other than issuing a warning/error
    #     (2) One problem seen is near-dups/remeasures that are still marked as 'P' in the database (use check-near-dups)
    #     (3) After that, think about doing some sort of grep of the flat-files
    for obs80_bit in db_only:
        R = process_db_only(cnx, obs80_bit, obs_db, DEBUG=DEBUG)
        if DEBUG and R.status != 1 :
            R.pprint()

    return {} 


def process_ff_only(cnx, obs80_bit, obs_ff, DEBUG=False ):
    ''' 
    # (b) For the N_ff_only observations, search the db based on obs80_bit to see whether these can be recovered
    #     Are they published?
    #     Are they listed under a different designation?
    #     What action/query would need to be taken to "fix" the database ?
    #     - Likely to not just be obs80 string, but also other fields
    '''
    # check input ok
    assert obs80_bit in obs_ff, 'bit not in obs dict (ff)'

    # object to hold results
    R = repair() 
    
    # Search db based on obs80_bit
    rows = psql.get_all_obs_by_obs80bit(obs80_bit, cnx)
    
    # If no data returned
    if rows is None :
        # Save 'status', 'action', 'info'
        R.app(-2, None, obs80_bit)
        
    # If there are returned rows
    else:
        for obs80_db, obsId, status, id in rows:

            # if status == 'X' and the entire obs80 matches ( or can be made to match )
            # =>  update status to 'P'/'p'
            if status in ['X','?'] and categorize_differences( obs_ff[obs80_bit], obs80_db , DEBUG = False).status == 0 :
                
                # Save 'status', 'action', 'info'
                # The *action* here should become some sql statement ...
                R.app( +1, f'Set status == P : {obs80_bit}', None)

            # some other category ...
            # elif ... :
            # pass 
 
            # anything remaining has not been explicitly dealt with 
            else :
                # Save 'status', 'action', 'info'
                R.app( -1, None, (obs80_db, obsId, status, id))
 
                    
    return R


def categorize_differences( obs_ff, obs_db , DEBUG = False) :
    ''' categorize the differences between two observations '''
    assert len(obs_ff) == len(obs_db), f'CANNOT PROCESS UNEQUAL OBS-LENGTHS: len(a)=len(obs_ff)={len(obs_ff)} , len)b_==len(obs_db)={len(obs_db)}'

    # Splitting up any 160 into 2x80 
    if len(obs_ff) > 80 :
        obs_ff_ = [obs_ff[:80],obs_ff[80:]]
        obs_db_ = [obs_db[:80],obs_db[80:]]
    else:
        obs_ff_ = [obs_ff]
        obs_db_ = [obs_db]

    # Categorizing differences we have considered & and are happy to deal with
    # These will only be done for "repair" tags that add info from ff to db: see below 
    # - NB Value tuple could be expanded to carry future 'operations' 
    add_dict = {
        (0,5)   : ('NUM') , # Number is missing
        (5,12)  : ('DES') , # Provisional Desig is missing
        (12,13) : ('AST') , # Asterisk
        (13,14) : ('PRO') , # Program Code
        (70,71) : ('BND') , # MagBand [I think]
        (71,72) : ('CAT') , # AstCat [I think]
        (72,73) : ('REF') , # Pubn Ref
        (72,74) : ('REF') , # Pubn Ref
        (72,75) : ('REF') , # Pubn Ref
        (72,76) : ('REF') , # Pubn Ref
        (72,77) : ('REF') , # Pubn Ref
        (74,77) : ('REF') , # Pubn Ref
        (71,77) : ('B+R') , # AstCat + PubRef
        (77,80) : ('COD') , # ObsCode
        (72,80) : ('R+C') , # PubRef + ObsCode
        (73,80) : ('R+C') , # PubRef + ObsCode
        (74,80) : ('R+C') , # PubRef + ObsCode
        (75,80) : ('R+C') , # PubRef + ObsCode
        (76,80) : ('R+C') , # PubRef + ObsCode
        }
    
    del_dict = {
        (70,71) : ('BND') , # Mag Band
    }
    
    # object to hold output
    R = repair() 

    # Loop over the 1 or 2 components of the 80/160 char string observation
    for obs_ff, obs_db in zip(obs_ff_, obs_db_):

        # If they are identical, don't need to analyze ...
        if obs_ff == obs_db:
            # Save 'status', 'action', 'info'
            R.app(0, None, None)

        else:
            # NB Sequence matcher tells you how to turn a --into--> b ...
            # ... so, as we want to make db --> ff, we should give db first ...
            opcodes = difflib.SequenceMatcher(None, obs_db, obs_ff).get_opcodes()
            tags    = np.array([_[0] for _ in opcodes])
            replace_indicees = np.where( tags == "replace" )[0]
            equal_indicees   = np.where( tags == "equal" )[0]
            other_indicees   = np.where( (tags != "equal") & (tags != "replace") )[0]
            assert len(replace_indicees) + len(equal_indicees) + len(other_indicees) == len(opcodes)

            print(obs_db[12:65] == obs_ff[12:65] , 
                  obs_db[65:] != obs_ff[65:] ,
                  obs_db[65:].count(' ') >= obs_ff[65:].count(' '),
                  obs_db[77:80] in ['C51','247'])
            # NB Have seen some oddities related to the way SequenceMatcher splits 70:80
            # ... (it can split it into odd strange sections when there are partial matches across PubRef & ObsCode, etc)
            # To handle a few different cases seen, I am allowing changes across 70-80
            #
            # *** THIS IS PRETTY AGGRESSIVE : I AM ALLOWING ANY REPLACEMENT IN THE LAST 10 CHARS ( AS LONG AS THERE ARE FEWER WHITE-SPACES) ***
            # *** *** SHOULD REVIEW *** *** 
            #
            if obs_db[12:70] == obs_ff[12:70] and \
               obs_db[70:] != obs_ff[70:] and \
               obs_db[70:].count(' ') >= obs_ff[70:].count(' ')    \
               or \
               obs_db[12:65] == obs_ff[12:65] and \
               obs_db[65:] != obs_ff[65:] and \
               obs_db[65:].count(' ') >= obs_ff[65:].count(' ') and  \
               obs_db[77:80] in ['C51','247']:

                
                # Save 'status', 'action', 'info'
                R.app(  +1 , obs_ff, None )
                
            # If we see any of these "other" indicees (outside of 70:80), immediately "give up" ...
            # ... because nothing yet coded for how to deal with these. 
            elif len(other_indicees) > 0 :
                
                # Save 'status', 'action', 'info'
                info_str = [ '{:7}   db[{}:{}] --> ff[{}:{}] \t:\t {!r:>8} --> {!r}'.format(
                    tag, i1, i2, j1, j2, obs_db[i1:i2], obs_ff[j1:j2]) for tag, i1, i2, j1, j2 in [opcodes[ind] for ind in other_indicees] ]
                info_str = f'CODE NOT YET DEVELOPED TO HANDLE THIS TAG' +'\t'*5 +'\n'+\
                    f'\tff:{obs_ff}\n'+\
                    f'\tdb:{obs_db}\n'+\
                    '\n'.join(info_str)
                R.app(  -2, None, info_str )  

                    

            # If all of the "repair" tags involve indicees in "add_dict", then I am happy to replace if ...
            # ... the same sections of the string are being dealt with in the db & ff strings [i1 == j1 and i2 == j2 ] ...
            # ... there is more/same information in the ff-section than the db-section [ obs_db[i1:i2].count(' ') >= obs_ff[j1:j2].count(' ')]
            # If these are satisfied, accept obs_ff as a replacement for obs_db
            elif  np.all(
                    np.array(
                        [ tag=='replace' and \
                          i1 == j1 and i2 == j2 and \
                          (i1,i2) in add_dict and  \
                          obs_db[i1:i2].count(' ') >= obs_ff[j1:j2].count(' ') \
                          for tag, i1, i2, j1, j2 in [opcodes[ind] for ind in replace_indicees] ])) :
                # Save 'status', 'action', 'info'
                # The *action* here should become some sql statement ...
                R.app(  +1 , obs_ff, None )
                    
 
            # Deletions from database
            # I am being deliberately cautious about this, so taking it slowly ...
            #
            # *** DO WE WANT TO ENCODE SOMETHING TO ALTER / UPGRADE THE FF-STRING ??? *** 
            # 
            elif  np.all(
                    np.array(
                        [ tag=='replace' and \
                          i1 == j1 and i2 == j2 and \
                          (i1,i2) in del_dict and  \
                          obs_db[i1:i2].count(' ') < obs_ff[j1:j2].count(' ') \
                          for tag, i1, i2, j1, j2 in [opcodes[ind] for ind in replace_indicees] ])) :
                # Save 'status', 'action', 'info'
                # The *action* here should become some sql statement ...
                R.app(  +1 , f'{obs_db} --> {obs_ff}', None )
            
                
            # Anything else is not yet encoded 
            else:
                
                # Save 'status', 'action', 'info'
                info_str = [ '{:7}   db[{}:{}] --> ff[{}:{}] \t:\t {!r:>8} --> {!r}'.format(
                    tag, i1, i2, j1, j2, obs_db[i1:i2], obs_ff[j1:j2]) for tag, i1, i2, j1, j2 in [opcodes[ind] for ind in replace_indicees] ]
                info_str = f'CODE NOT YET DEVELOPED TO HANDLE THIS DIFFERENCE' +'\t'*5 +'\n'+\
                    f'\tff:{obs_ff}\n'+\
                    f'\tdb:{obs_db}\n'+\
                    '\n'.join(info_str)

                R.app(  -1, None, info_str )
                                      
    
    # return after combining / summarizing the necessary changes / categorization
    #print(R.status, R.action, R.info)
    R.summarize()
    return R


def process_db_only(cnx, obs80_bit, obs_db, timeDeltaSeconds=20, arcsecRadius=5.0, DEBUG=False ):
    ''' 
    # (d) For the N_db_only observations: what the hell do we do about these? 
    #     (0) Start by doing nothing other than issuing a warning/error
    #     (1) One problem seen is near-dups/remeasures that are still marked as 'P' in the database (use check-near-dups)
    #     (2) After that, think about doing some sort of grep of the flat-files
    '''

    # check input ok
    assert obs80_bit in obs_db, 'bit not in obs dict (db)'

    # object to hold results
    R = repair() 

    # (1) Check whether this problematic obs80_bit is near-dup of any other obs80_bit from db
    nds = cnd(obs80_bit, obs_db, timeDeltaSeconds=timeDeltaSeconds, arcsecRadius=arcsecRadius, DEBUG=DEBUG )
    if len(nds) :
        # Save 'status', 'action', 'info'
        R.app(  +1 , f'NEAR-DUP TO REMOVE:{obs80_bit}\t\t\t\nDUP:{obs_db[obs80_bit]} OK_:{obs_db[nds[0]]}', None )
        
    # (2) ... other ????
    elif False :
        pass
    
    # (X) *Not fixed* (no diagnosis / fix coded-up) 
    else:
        # Save 'status', 'action', 'info'
        R.app(  -1 , None, obs80_bit)
        
    return R 

def cnd(obs80_bit, obs_db, timeDeltaSeconds=20., arcsecRadius=5.0, DEBUG=False):
    ''' check near dups '''
    # Ddoing a version of check-near-dups (but without having to query db again)
    # Because db-only is rare, I will do obs80 conversion to facilitate comparison
    #  - N.B. This is a hack-approach, because I have double-length lines here ...
    # Select any observations that are close in time
    #  - The allowed time-range is the time of the observation +/- the time-delta (in days)
    delta_days = timeDeltaSeconds / (86400.)
    try:
        target     = Opt(obs_db[obs80_bit])
        obs80_dict = { k: Opt(line) for k, line in zip(list(obs_db.keys()),obs_db.values()  ) if np.abs( Opt(line).jdutc - target.jdutc) <  delta_days and k != obs80_bit}
    except:
        print('PROBLEM PARSING w/ Opt ')
        for k,v in obs_db.items():
            
            try:
                Opt(obs_db[k])
            except:
                print(k,':::',v)
        print('PROBLEM PARSING w/ Opt ')
        sys.exit()

    # Filter results to only those within 'arcsecRadius' of the supplied originalObs80-string
    near_dups = []
    if obs80_dict:
        dotProducts         = np.clip([ hp.ang2vec(np.radians(90 - target.dec), np.radians(target.ra*15.)).dot(uv) for uv in \
                                      [ hp.ang2vec(np.radians(90 - o.dec)     , np.radians(o.ra*15.)) for o in obs80_dict.values() ] ], -1,1)
        angles              = np.degrees(np.arccos(dotProducts))*3600.
        #indicees           = np.where( angles < arcsecRadius )[0]
        near_dups           = np.array(list(obs80_dict.keys()))[ np.where( angles < arcsecRadius )[0] ]
    
    return near_dups
        
def Opt(line):
    ''' Hacked version of sonia's Optical object (always taking first 80 chars ...) '''
    return o80.Optical(
        num=line[:5].strip(),
        desig=line[5:12].strip(),
        disc=line[12].strip(),
        note1=line[13].strip(),
        note2=line[14].strip(),
        jdutc=o80.date2JD(line[15:32]),
        ra=o80.RA2hrRA(line[32:44]),
        dec=o80.Dec2degDec(line[44:56]),
        mag=o80.floatMag(line[65:70]),
        band=line[70].strip(),
        cod=line[77:80]
    )

    
if __name__ == "__main__":
    #check_consistency(int(sys.argv[1]),
    #                  int(sys.argv[2]),
    #                  DEBUG = True if len(sys.argv) > 3 and ( bool(sys.argv[3]) or sys.argv[3] == 'DEBUG') else False)
    
    
    # (i) Start by checking the internal self-consistency of the flat-files
    check_flat_file_internal_consistency(   int(sys.argv[1]),
                                            int(sys.argv[2]),
                                            DEBUG = True if len(sys.argv) > 3 and ( bool(sys.argv[3]) or sys.argv[3] == 'DEBUG') else False )
    
