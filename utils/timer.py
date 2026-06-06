from time import time

class Timer:
    timers = {}
    @staticmethod
    def attach_timer(name="timer"):
        Timer.timers[name] = 0
        def decorate(func):
            def func_wrapper(*args, **kwargs):
                start_time = time()
                res = func(*args, **kwargs)
                Timer.timers[name] +=  time() - start_time
                return res
            return func_wrapper
        return decorate
    
    