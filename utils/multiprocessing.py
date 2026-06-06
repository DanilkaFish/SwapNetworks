from __future__ import annotations

from typing import Callable


import multiprocessing as mp

class Threads:
    def __init__(self):
        self.manager = mp.Manager()
        # data = manager.dict()
        self.procs = []
        pass

    def get_shared_dict(self):
        return self.manager.dict()

    def set_parallel_function(self, fun: callable):
        self.fun = fun

    def add_proccess(self, args):
        self.procs.append(mp.Process(target=self.fun, args=args))

    def start_join_close(self):
        self.start()
        self.join()
        self.close()

    def start(self):
        for proc in self.procs:
            proc.start()

    def join(self):
        for proc in self.procs:
            proc.join()

    def close(self):
        for proc in self.procs:
            proc.close()