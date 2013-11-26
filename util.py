import time
from functools import wraps


def time_it(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        t0 = time.time()
        out = f(*args, **kwargs)
        elapsed_time = time.time() - t0
        print '-- %s took %s sec to run' % (str(f.__name__), str(elapsed_time))
        return out
    return wrapper
