    #!/usr/bin/env python3
'''

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
from collections import defaultdict, Counter, OrderedDict
import textwrap 
import healpy as hp
import glob
import re
import shutil
import copy
from functools import lru_cache

# ------------------------------------------------------------------------
# RAY/DASK PARALLELIZATION
# ------------------------------------------------------------------------
#import ray
import dask
from distributed import Client

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
# HIGH-LEVEL CONSISTENCY CHECKING FUNCTIONS
# ------------------------------------------------------------------------
@dask.delayed
def establish_internal_consistency_of_flat_files_for_single_desig( desig, proc_dir, cnx=None, DEBUG = True ):
    '''
    '''
    print(f'Inside *establish_internal_consistency_of_flat_files_for_single_desig*: desig={desig}')
    
    # If we did not get a connection passed in, establish one
    # This is particularly important for ray-parallelization
    if cnx is None :
    
        # Establish one-off connection to the database on mpcdb1
        # NB Despite being flat-file focused, we might need to
        # get submission IDs from the database ...
        host     = 'mpcdb1.cfa.harvard.edu'
        database = 'vmsops'
        cnx = psycopg2.connect(f"host={host} dbname={database} user=postgres")

    
    # Get obs from flat files
    # - returns a list of obs
    print('Getting obs from ff for ', desig)
    obs_ff = ff.get_obs_from_ff(desig, proc_dir, DEBUG=DEBUG)

    '''
    # Look for duplicates
    deduped_obs_list, duplicates = ff.find_duplicates(obs_ff)

    # Look for 2-line obs
    # Combine together where possible
    # Identify problems where exist
    obs_dict, orphans = ff.combine_two_line_obs(deduped_obs_list)

    # Are there other problems with the flat-file data that we can look for ?
    # (1) - Sometimes we do not have "Note 2" before 2020 in obs80: replace blank with default ?

    # Identify and fix any problems within the flat-file data
    if duplicates or orphans :
        
        # Fix simple duplicates ...
        if duplicates:
            print('Fixing duplicates ...', desig)
            report = ff.fix_primary_flat_file_data(desig, duplicates, [] , DELETING=True )
           
        # Attempt to fix sat/roving stuff here ...
        if orphans:
            print('Fixing orphans ...', desig)
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
            #if incorrect_list != []:
            #    report = ff.fix_primary_flat_file_data(desig, incorrect_list, correct_list )
            #    print(f' report from fix_primary_flat_file_data : {report} ')
    else:
        print(f'\t ... no problems found for {desig}')

    
    '''



# ------------------------------------------------------------------------
# FUNCTIONS RELATED TO READING FROM FLAT-FILES
# ------------------------------------------------------------------------
def get_obs_from_ff(desig, proc_dir, DEBUG=False):
    '''
    # Get obs from flat files 
    # - Do we need to search using findn & findu on multiple "identifications" ?  
    # - Would be good to double-check nothing extra is returned 
    # - What files do / don't get searched?
    #
    #
    # Is findn automatically sorted ?
    #
    '''
    if DEBUG:
        pass # print('ff-query is basic: needs to be expanded to use findu & other (linked) designations')
        
    # Query stolen from Dave (see 'getobs_num.py)
    obs = []
    outfile = os.path.join(proc_dir, desig + "_ff.obs")
    cmd = f"""$findn -o={outfile} -O={desig}"""
    os.system(cmd)
    if os.path.isfile(outfile):       
        with open(outfile) as f:
            lines = f.readlines()
            obs_list = [_.strip() for _ in lines]
    else:
        print(f'outfile:{outfile} did *NOT* get generated')
            
    # Explicitly sort by time ?
    
    # Get rid of the file
    # This is to stop the build-up of huge quantities of data in
    # the local and/or processing directories.
    os.remove(outfile)
    
    return obs_list

def find_duplicates(obs_list):
    ''' at the moment this is targeted at obs for a single object...
        but there's nothing to stop it being used on obs from across many objects ...
        as it uses obs80-bit
    '''
    # Find any duplicate lines
    obs = OrderedDict()
    for n,line in enumerate(obs_list):
        obs80_bit       = line[15:56]
        if obs80_bit in obs:
            obs[obs80_bit].append(line)
        else:
            obs[obs80_bit] = [line]
    
    # Decide which line of any duplicates will be kept
    deduped_obs_list, probs = [], []
    for obs80_bit, lineList in obs.items():
        # Singles are fine
        if len(lineList) ==1:
            deduped_obs_list.extend(lineList)
        
        # Put duplicates into problem-list for now
        # I'm choosing one to keep by sorting on the pubn-ref and keeping the earliest
        # If they are the same, it doesn't matter as they are duplicates
        else:
            pub_refs = [_[72:77] for _ in lineList]
            pub_refs, lineList = zip(*sorted(zip(pub_refs, lineList)))
            
            # Keep one, discard the rest
            deduped_obs_list.append(lineList[0])
            probs.extend( lineList[1:] )

    return deduped_obs_list, probs
    
    
def combine_two_line_obs(deduped_obs_list):
    '''
    # Find and combine 2-line obs ...
    # ... make dict key-ed on obs80 bit
    # NB : The check is being made under the
    # assumption that any dups have been removed ...
    '''
    obs_dict, probs = {}, []
    for n,line in enumerate(deduped_obs_list):
        obs80_bit       = line[15:56]

        # Satellites & Roving
        if line[14] in ['S','V']   and deduped_obs_list[n+1][14] == line[14].lower()  :
            obs_dict[obs80_bit] = line + deduped_obs_list[n+1]
        elif line[14] in ['s','v'] and deduped_obs_list[n-1][14] == line[14].upper()  :
            pass # because of previous logic
            
        # Could split these problems if desired / useful
        elif line[14] in ['S','V'] and deduped_obs_list[n+1][14] != line[14].lower()  or \
             line[14] in ['s','v'] and deduped_obs_list[n-1][14] != line[14].upper()  :
            probs.append(line)

        # Standard lines 
        else:
            obs_dict[obs80_bit] = line
            
    return obs_dict, probs
    

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

    NB1 : I wanted to use formatobs.txt, but sometimes doesn't exist, so using "artifact" instead
    NB2 : 3317 does not have a "formatobs" file ... gaps in artifact file ... how did it get processed ?

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

    The assumption is that
    (i) the "incorrect_published_obs80" string
    has been used to find the original submission-artifact
    (via *find_original_submission_artifact()* ),
    and
    (ii) the original_obs80 extracted
    (using *extract_original_observation()* ).

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

    # If initial version of correct_obs80 doesn't exist, then crap out as can't fix ...
    try:
        print('extended___________________ ', correct_obs80)
    except:
        sys.exit('cannot fix ...', original_obs80 , incorrect_published_obs80)

    # Using the cat_diff routine that was original developed to compare flat-files / db
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
    *** (ii) Loop over any files that need to be fixed ...
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
    
    #*** (i) Find the relevant primary data file [can be in /sa/mpn or in tot*]
    src_files = []
    for incorrect_published_obs80 in incorrect_list:
        src_files.extend( find_primary_data_file( desig, incorrect_published_obs80 ) )
        src_files = list(set(src_files))
    print('src_files = ', src_files)
    assert len(src_files), f'No src_file could be found that contains the incorrect data ... incorrect_list={incorrect_list}'
    
    #*** (ii) Loop over any files that need to be fixed ...
    dst_dir = newsub.generate_subdirectory( "obs_cons" )
    for src_file in src_files:
        fix_single_file(src_file, dst_dir, desig, incorrect_list, correct_list, DELETING=DELETING)
        
    
    return True
        
def fix_single_file(src_file, dst_dir, desig, incorrect_list, correct_list, DELETING=False):
    '''

    *** Need to be really careful about this ***
    *** Need to do something like ... ***
    *** (i) Freeze the system (lock status)
    *** (ii) Copy primary data file to temp location (esp. while developing)
    *** (iii) Find the location of the incorrect data in the primary data file
    *** (iv) Replace the incorrect data with the correct data (non-trivial : needs to have pubn-record, etc)
    *** (v) Do some sense checks of the difference between the initial and fixed versions
    *** (vi)  write the data to the temp file
    *** (vii) Replace the primary data with the fixed copy
    *** (viii) Unlock the system

    inputs:
    -------

    returns:
    --------

    '''
    # We want to 'permanently' save some output files ...
    save_dir = '/sa/conchecks/data_products/'


    # (i) Freeze the system (lock status)
    # ~~~~~~~~~~~~~ IF WE CRAP OUT AT ANY POINT BELOW WE NEED TO RELEASE THE LOCK ~~~~~~~~~~~~~~~~~~
    print(f'{desig}: Checking & Setting mpc_temp_status')
    desired_status_string = f'Fixing Flat Files : {desig} : {src_file}'
    while mpc_status.get_status("mpc_temp_status") != desired_status_string :
        print(f'{desig} waiting for mpc_temp_status to be released ...')
        time.sleep(np.random.rand()*0.01)
        if mpc_status.get_status("mpc_temp_status") == '':
            mpc_status.set_status("mpc_temp_status", desired_status_string)
        
    # If we get here, then it must have been set to desired_status_string,
    # so this process has the lock and can continue
    assert mpc_status.get_status("mpc_temp_status") == desired_status_string
    print('I now have the status lock ... ', desired_status_string )
    
    
    # Wrapping everything in try-except loop with the goal of ensuring the status-lock is always released.
    try:
    
        #*** (ii) Copy primary data file to temp location
        #(esp. important while developing)
        dst_file= os.path.join(  dst_dir , os.path.basename(src_file) ); print('dst_file = ',dst_file)
        shutil.copyfile(src_file, dst_file )

        # Read the primary data file
        with open(dst_file,'r') as fh : data = fh.readlines()

        # Files to write to so that MR can update mysql
        bad_filepath  = os.path.join(save_dir, desig + '_bad.dat')
        good_filepath = os.path.join(save_dir, desig + '_good.dat')
        print(f' Files for MYSWL: bad_filepath= {bad_filepath} , good_filepath= {good_filepath} ')
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
                        assert isinstance(incorrect_published_obs80, str), \
                            f'incorrect_published_obs80 is of type {type(incorrect_published_obs80)}, rather than a string'

                        #*** (iv) Find the location of the incorrect data in the primary data file
                        line_num = [i for i,line in enumerate(data) if incorrect_published_obs80.strip() in line]
                        assert len(line_num) < 3, \
                            'len(line_num)={len(line_num)} which is >=3 which seems like a suspiciously large number so I am terminating...'

                        #*** (v) Replace the incorrect data with the correct data (the correct data has been created earlier)
                        #        At the same time we also output the incorrect & correct data to some files to be used to update the MYSQL database
                        for n,line in enumerate(data):
                            if n not in line_num :
                                # We keep the normal stuff as-is
                                fixed_data.append(line)
                            else:
                                # For removal from mysql
                                bad_fh.write(line)

                                # 1-line ...
                                if isinstance(corrected_obs80, str):
                                    l = corrected_obs80 if corrected_obs80[-1]=='\n' else corrected_obs80+'\n'
                                    # Corrected data for flat files
                                    fixed_data.append(l)
                                    # Corrected data for mysql
                                    good_fh.write(l)
                                # 2-line ...
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
                        assert len(fixed_data) - len(data) == len(line_num), \
                            'Lengths do not make sense: {len(fixed_data),len(data),len(line_num)} '

                        # copy fixed data into data ready for next loop around ...
                        data = copy.deepcopy(fixed_data)
                    

        #*** (vii)  write the data to the temp file
        #           NB This deals with both the DELETIONS and the FIXES
        replace_file = dst_file + 'replace'
        assert not os.path.isfile(replace_file), \
            'replacement file {replace_file} already exists which is bad'
        with open(replace_file,'w') as fh :
            for line in fixed_data:
                l = line if line[-1]=='\n' else line+'\n'
                fh.write(l)
        assert os.path.isfile(replace_file), \
            'replacement file {replace_file} does NOT exist which is bad'


        #*** (viii) Replace the primary data with the fixed copy
        print('THE FOLLOWING COMMAND HAS BEEN COMMENTED-OUT !!!')
        print(f'replacing file={src_file} with file {replace_file} ')
        #shutil.copyfile(replace_file, src_file)
        print('THE ABOVE WAS NOT DONE!!')
        
        
        #*** (ix) Recreate the index files if necessary
        #         NEED TO BE CAREFUL ABOUT THIS ...
        #         (a) Mike/Dave indicated this is only necessary if the file being altered is one of the permanent,
        #             master files, rather than one of the temp *tot* files
        #         (b) However, my inspection of /share/apps/mpec/publish_dou_mpec.sh, /share/apps/com/indexed/update.sh
        #             (and sub-scripts) suggests that there *ARE* some form of index files for the temp/pending/within-month files
        #         (c) TO gain some understanding, the monthly-prep rebuilds are done here : /sa/com/indexed/update.sh [SAME AS ABOVE]
        #         (d) Given that ... calls /share/apps/com/indexed/buildnumupd.sh
        #   *** I recall reading and re-reading the scripts that remake the numbered index files ***
        #   *** And I recall deciding that perhaps it wasn't necessary to IMMEDIATELY recreate the index files ***
        #   *** This seems counter-intuitive : but the reasoning is as follows ...  ***
        #   *** The index files do NOT actually point to the primary data files, but instead point to a COPY (without line-breaks) ***
        #   *** This means that changes to the primary data do NOT cause malfunctions in the indexing ***
        #   *** Hence we can do a whole raft of changes to the primary data files, and then only do a single re-indexing ***
        #   *** This is useful as the reindexing takes a long time ***
        #   *** NB, this obviously leaves the index files temporarily serving up crappy data ... ***
        #   ***     but I regard this as acceptable as they were serving up crappy data before any fix anyway ***

        #*** (xi) Remove / Tidy-up the temp files & temp dir
        #shutil.rmtree(dst_dir)
        #assert not os.path.isdir(dst_dir), f'dst_dir={dst_dir} still exists when it should NOT'

        #*** (ix) Unlock the system
        print('Unsetting the mpc_temp_status')
        mpc_status.set_status("mpc_temp_status","")

    
    except Exception as e: 
        print('\n'*2)
        print('EXCEPTION IN fix_single_file')
        print('\n'*2)
        print(e)
        print('\n'*2)
        print('Unsetting the mpc_temp_status as part of the EXCEPTION handling')
        if mpc_status.get_status("mpc_temp_status") == desired_status_string:
            mpc_status.set_status("mpc_temp_status","")

    return True

@lru_cache(16)
def get_potential_source_files(desig):
    ''' Get filenames of potential source files for numbered object '''
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
    
    return potential_source_files
    
@lru_cache(16)
def get_contents_of_potential_source_files(desig):
    ''' read the contents of potential source files for numbered object '''

    # Get names of potential source files
    potential_source_files = get_potential_source_files(desig)

    # Read contents into a dict
    content_dict = {}
    for f in potential_source_files:
        with open(f, "r") as fh:
            content_dict[f] = fh.readlines()

    return content_dict
    
def find_primary_data_file( desig, incorrect_published_obs80 ):
    ''' #*** (i) Find the relevant primary data file [can be in /sa/mpn or in tot*] '''
    
    # Get contents of potential source files
    content_dict = get_contents_of_potential_source_files(desig)
    
    # Check which file(s) incorrect_published_obs80 is in
    src_files = []
    for f, data  in content_dict.items():
        for line in data:
            if incorrect_published_obs80.strip() in line:
                src_files.append(f)
                    
    # NB I am deliberately *NOT* returning the line number here ...
    # ... because there is a TINY chance the file will be altered ...
    # ... between now and when I lock the file 
    return list(set(src_files))
   

# ------------------------------------------------------------------------
# RAY PARALLELIZATION (FOR CROSS-DESIGNATION DUPLICATION CHECKS)
# ------------------------------------------------------------------------
@dask.delayed
def check_two_files_for_dups(f1,f2, i1,i2):
    ''' '''
    print(f'check_two_files_for_dups:{f1,f2, i1,i2}')
    
    # Read first file into a dict
    with open(f1,'r') as fh1:
        d1 = {line[15:56]:True for line in fh1 if line[14] not in ['s','v']}

    # Go through the lines in the second file and for any duplicates, ...
    # ... add *both* integers to a dups dict
    DUPS = {}
    with open(f2,'r') as fh2:
        for line in fh2:
            key = line[15:56]
            if key in d1 :
                DUPS[key]=[i1,i2]
    del d1
    return DUPS
    
@dask.delayed
def compare_file_against_provided_file_dict(d1,f2, i1,i2):
    ''' '''
    print(f'compare_file_against_provided_file_dict:{f2, i1,i2}')

    # Go through the lines in the second file and for any duplicates, ...
    # ... add *both* integers to a dups dict
    DUPS = {}
    with open(f2,'r') as fh2:
        for line in fh2:
            key = line[15:56]
            if key in d1 :
                DUPS[key]=[i1,i2]
    del d1
    return DUPS

#@dask.delayed
def compare_two_dicts_for_dups(d1,d2, i1,i2):
    ''' '''
    print(f'compare_two_dicts_for_dups:{i1,i2}')

    # intersecn indicates duplicate obs80-bits
    intersecn = d1.keys() & d2.keys()

    return {k:[i1,i2] for k in intersecn}