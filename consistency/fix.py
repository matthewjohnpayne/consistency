'''
MJP
2020-12-13

Use this to fix any problems identified by the monitoring script.

'''


# imports 
import sys, os 
import flat_file_duplicates as ffd
 
# fix duplicates across flat_files ...
ffd.CrossDesignationDuplicates.fix() 


