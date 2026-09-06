import os
from time import time as time_now, ctime

def non_negative_args(fn):
    def check_non_negative(*args, **kwargs):
        for a in args:
            if isinstance(a, (int, float)):
                if a < 0:
                    try:
                        with open('error_log.csv', 'x') as f:
                            f.write("user,timestamp,message\n")
                            f.close()
                    except:
                        pass
                    finally:
                        with open('error_log.csv', 'a') as f:
                            errmsg = "Error: Input argument must be non-negative.\n"
                            error_entry = os.getlogin() + ',' + ctime(time_now()) + ',' + errmsg
                            f.write(error_entry)
                            f.close()      
                        raise ValueError(errmsg)
                        
        for k, kv in zip(kwargs.keys(), kwargs.values()):
            if isinstance(kv, (int, float)):
                if kv < 0:
                    try:
                        with open('error_log.csv', 'x') as f:
                            f.write("user,timestamp,message\n")
                            f.close()
                    except:
                        pass
                    finally:
                        with open('error_log.csv', 'a') as f:
                            errmsg = "Error: Input {} must be non-negative.\n".format(k)
                            error_entry = os.getlogin() + ',' + ctime(time_now()) + ',' + errmsg
                            f.write(error_entry)
                            f.close()      
                        raise ValueError(errmsg)
        return fn(*args, **kwargs)
    return check_non_negative
