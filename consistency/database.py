#!/usr/bin/env python3
'''

'''

# ------------------------------------------------------------------------
# THIRD PARTY IMPORTS
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# MPC IMPORTS
# ------------------------------------------------------------------------
import mpc_psql as psql


    
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
