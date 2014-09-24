# Built-in modules #
from datetime import datetime

################################################################################
class Timer(object):
    """Useful for timing the different steps in a pipeline"""

    def __init__(self):
        self.start_time = datetime.now()
        self.last_mark = datetime.now()

    def print_start(self):
        print "Start time: %s" % self.start_time

    def print_elapsed(self, reset=True):
        print "Elapsed time: %s" % (datetime.now() - self.last_mark)
        self.last_mark = datetime.now()