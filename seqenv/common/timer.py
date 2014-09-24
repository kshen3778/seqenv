# Built-in modules #
from datetime import datetime

# Internal modules #
from seqenv.common.color import Color

################################################################################
class Timer(object):
    """Useful for timing the different steps in a pipeline"""

    def __init__(self):
        self.start_time = datetime.now()
        self.last_mark = datetime.now()

    def print_start(self):
        print self.prefix + "Start time: %s" % (self.start_time) + self.suffix

    def print_elapsed(self, reset=True):
        print self.prefix + "Elapsed time: %s" % (datetime.now() - self.last_mark) + self.suffix
        self.last_mark = datetime.now()

    @property
    def color(self):
        """Should be use color or not ? If we are not in a shell, then not"""
        import __main__ as main
        if not hasattr(main, '__file__'): return True
        return False

    @property
    def prefix(self): return "" if not self.color else Color.grn
    @property
    def suffix(self): return "" if not self.color else Color.end
