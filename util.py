import time
from threading import Thread
from multiprocessing import Process
from functools import wraps

from Queue import Queue

def time_it(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        t0 = time.time()
        out = f(*args, **kwargs)
        elapsed_time = time.time() - t0
        print '-- %s took %s sec to run' % (str(f.__name__), str(elapsed_time))
        return out
    return wrapper


class Pool():

    def __init__(self, num_threads):
        self.queue = Queue()
        self.num_threads = num_threads
        #for i in range(num_threads):
        #    t = Process(target=self._run_thread)
        #    t.daemon = False
        #    t.start()

    def apply_async(self, f, args):
        q = Queue()
        self.queue.put((q, f, args))
        return q

    def apply_async1(self, f, args):
        q = Queue()
        p = Process(target=f, args=args)
        p.start()
        class A():
            def get(self):
                p.join()
                return
        return A()

    def _run_thread(self):
        print 'thread waiting'
        while True:
            q, f, args = self.queue.get()
            q.put(f(*args))


    def concurrent(self, f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            q = Queue()
            self.queue.put((q, f, args))
            return q.get()
        return wrapper
