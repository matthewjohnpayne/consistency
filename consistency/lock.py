# ------------------------------------------------------------------------
# THIRD PARTY IMPORTS
# ------------------------------------------------------------------------
import ray
from filelock import FileLock

# ------------------------------------------------------------------------
# MPC IMPORTS
# ------------------------------------------------------------------------
import status as mpc_status

# ------------------------------------------------------------------------
# Class/Actore to facilitate locking / acquiring "mpc_temp_status"
# ------------------------------------------------------------------------
@ray.remote
class Locker:
    timeout = 1000
    
    def __init__(self):
        # Using a separate lock file to prevent my multiple
        # workers/functions stealing it from one another
        self.lock = filelock.FileLock(os.path.expanduser("~/.sim.lock"))

    def acquire_status(self, desired_status_string):
        # Don't need to do anything if we already have the status !!
        if mpc_status.get_status("mpc_temp_status") != desired_status_string:
        
            try:
            
                # wait to acquire lock from parallel workers
                with self.lock.acquire(timeout=timeout):
                
                    # wait to aquire lock from any other code
                    # e.g. PV uses this during identification/linking
                    while mpc_status.get_status("mpc_temp_status") != desired_status_string:
                    
                        time.sleep(np.random.rand()*0.01)
                        
                        if mpc_status.get_status("mpc_temp_status") == '':
                            mpc_status.set_status("mpc_temp_status", desired_status_string)

                        time.sleep(np.random.rand()*0.01)
                        
                    assert mpc_status.get_status("mpc_temp_status") == desired_status_string, \
                        f'Problem: mpc_temp_status = {mpc_status.get_status("mpc_temp_status")}'

            except Exception as e:
                print('Problem with *aquire_status()*')
                print(e)
                print('\t:',desired_status_string)
        
        return mpc_status.get_status("mpc_temp_status")
        
    def relinquish_status(self, desired_status_string):
        
        # Only attempt to change if the temp_status is what you think it is ...
        if mpc_status.get_status("mpc_temp_status") == desired_status_string:
        
            # wait to acquire lock from parallel workers
            with self.lock.acquire(timeout=timeout):
            
                # set empty status
                mpc_status.set_status("mpc_temp_status", "")
                return True
        else:
            return False 
        
