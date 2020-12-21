'''
MJP 2020-12-13

Intend to use this as a *control* file for monitoring for problems
 - We'll want to monitor flat-files & obs-table, looking for a variety of different problems

I'll probably define the functions elsewhere and import them into here. 
 - A variety of WIP functions have been written in a few adjacent files, but I want to restructure/tidy things. 

'''

# imports 
import sys, os 
import flat_file_duplicates as ffd
 
# monitor for duplicates across flat_files ...
ffd.CrossDesignationDuplicates.find() 



