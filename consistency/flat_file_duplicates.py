#!/usr/bin/env python3
'''
MJP 2020-12-13
Some code to contribute towards establishing observational consistency
Intend to use this code purely for finding duplicates in flat-files
This will have to involve 2 types of duplicate
(i) duplicates within files
 - Not yet implemented within this file
(ii) duplicates across files
 - See "CrossDesignationDuplicates" class

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
from collections import Mapping, Container
from sys import getsizeof


# ------------------------------------------------------------------------
# RAY/DASK PARALLELIZATION
# - MJP 2020-12-13 : Not using these at the moment : just going to keep it slow and simple for now
# ------------------------------------------------------------------------
#import ray
#ray.init('auto')
##import dask
##from distributed import Client
##import dask.bag as db
##client = Client('tcp://131.142.192.121:8786')

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


# ------------------------------------------------------------------------
# CROSS DESIGNATION DUPLICATES...
# ------------------------------------------------------------------------
    
class CrossDesignationDuplicates():
    '''
    There's a possibility that the same observation has
    been published against multiple object-designations.
     - I.e. the same observation may appear in multiple flat-files
     - That is bad
    Let's check for duplicate observations appearing in different files
     - Numbered & Unnumbered
    
    NB I am not explicitly checking for duplicates WITHIN files here:
     - I am assuming I do that elsewhere ...
    '''

    def __init__(self, ):
    
        # We want to 'permanently' save some output files ...
        save_dir = newsub.generate_subdirectory( 'obs_cons' )
        print(f'CrossDesignationDuplicates saving into {save_dir}')
        
        self.map_file = filepath=os.path.join(save_dir , 'mapping.txt')
        self.dup_file = filepath=os.path.join(save_dir , 'duplicates.txt')


    def find(self, METHOD='ALL'):
        '''
        Find any duplicates
        
        There's more than one way to do this.
        The fundamental problem relates to having to check EVERY observation against EVERY OTHER ONE
        In this *find* function I give the option of calling two different methods ...
        '''
        
        # Get the files that need to be searched through
        file_dict = self._get_filenames()
        
        # ------- SEARCH FOR DUPLICATES -----------
        
        # (1) In this approach we read all files into memory at once
        if METHOD == 'ALL':
            duplicate_dict = self._read_all(file_dict)
        
        # (2) Alternatively, try to split the files up so that I don't have to read them all at once
        # - Nothing working yet
        elif METHOD == 'SPLIT':
            pass
            #duplicate_dict = self._read_split()
            
        # (3) Try some parallel read method
        # - Nothing working yet
        elif METHOD == 'PARALLEL':
            pass
            # self._find_duplicates_using_parallel_method()
        
        # (4) Explicitly fail if METHOD not in ['ALL','SPLIT','PARALLEL']
        else:
            print(f'METHOD not understood: METHOD = {METHOD}')
            
        # ------- SAVE ANY DUPLICATES TO FILE -----------
        self.save(duplicate_dict)
        
        
        
    def save(self, duplicate_dict):
        ''' save any duplicates to file '''
        with open( self.dup_file, 'w') as fh:
            for obs80bit, lst in DUP.items():
                for i,n in enumerate(lst):
                    fh.write(f'{obs80bit},{i},{file_dict[n]},{num[n]}\n')
        print('\t'*3,'created/updated:', self.dup_file)


        
    def _get_filenames(self, ):
        ''' get a dict containing all the filenames we want to work with ...'''

        # ------------ NUMBERED FILES ------------------
        # Primary, published files
        files_ = glob.glob(f'/sa/mpn/N*dat', recursive=True)
        
        # In-progress ( between monthly pubs) files are in different location ...
        # E.g. "tot.num", "pending.num", ..., ...
        files_.extend( glob.glob(f'/sa/obs/*num', recursive=True) )
        
        
        # ---------------- UN-numbered FILES -----------
        files_.extend(glob.glob(f'/sa/mpu/*dat', recursive=True))
        
        # *** WHY IS THERE NO "/sa/obs/*unn"??? ***
        #files_.extend( glob.glob(f'/sa/obs/*unn', recursive=True) )


        # ---------------- File-Mapping -----------
        file_dict = { n:f for n,f in enumerate(files_)}
        num       = { n:True for n,f in file_dict.items() } # Later on might want unnum files as well
        with open( self.map_file,'w') as fh:
            for n,f in file_dict.items():
                fh.write(f'{n},{f},{num[n]}\n')
        print('created...', self.map_file )
        
        return file_dict
    

    def _read_all(self, file_dict):
        '''---------------- Big data read ----------
        Read the data into a single, massive dictionary
        This is going to be challenging
        '''
        ALL = {}
        DUP = defaultdict(list)
        
        # Loop through all of the files ...
        for n,f in file_dict.items():
            print(n,f, len(file_dict))
            
            # Read the file contents into a "local" dictionary
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
            local     = {line[15:56]:n for line in fh if line[14] not in ['s','v']}
            
            del local
            del intersecn
            
            # Because I am impatient, I will print out the entire dict any time there is content ...
            if DUP:
                self.save(duplicate_dict)
            
        del ALL
        return duplicate_dict
    """
    def _read_split(self, file_dict):
        '''
        Not yet implemented
        The idea is that *_read_all* is slow, and that part of the problem is driven by having to read all files into memory sequentially.
        Perhaps we could implement an alternative method that splits the files up into "sets of files"
        And that then uses a "trivial parallelization over multiple machines" to read the files.
        This will end up using more processor-time (reading each file multiple times)
        But it would use less wall time (because of parallelization over machines)
        '''
        pass
    """
    

    """
    def _find_duplicates_using_parallel_method(self,file_dict):
        '''
        The idea is to use a parallel-engine like "ray" or "dask" to get the file-reads going in parallel
        I.e. reading across 16 (or 300) cores
        Commented-out as not yet fully developed/tested
        '''

        # ---------------- Big data read ----------
        print('reading ...')
        start = time.time()

        @dask.delayed
        def read_into_bit_dict(f):
            with open(f,'r') as fh:
                return {line[15:56]:True for line in fh if line[14] not in ['s','v']}


        list_of_dicts = [ read_into_bit_dict(f) for f in list(file_dict.values())[:100] ]
        #db.from_sequence(list(file_dict.values())[:10]).map(read_into_bit_dict)
        list_of_dicts = dask.compute(*list_of_dicts)
        print('...', len(list_of_dicts), type(list_of_dicts), [type(_) for _ in list_of_dicts], time.time()-start )
        sys.exit()
        list_of_dup_dicts = []
        for i,di in enumerate(list_of_dicts):
            for j, dj in enumerate(list_of_dicts[i+1:]):
                list_of_dup_dicts.append( ff.compare_two_dicts_for_dups(di,dj, i,j) )
                print('\t',len(list_of_dup_dicts[-1]))
        sys.exit()
    """

    def fix(self, identifier ):
        ''' Fix any duplicates that were previously found in the *find* routine (above)'''
        pass
