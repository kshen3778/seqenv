# Built-in modules #
import tempfile, datetime, dateutil

# One liners #
flatter = lambda x: [item for sublist in x for item in sublist]

################################################################################
def new_temp_path(**kwargs):
    """A new temporary path"""
    handle = tempfile.NamedTemporaryFile(**kwargs)
    path = handle.name
    handle.close()
    return path

################################################################################
class GenWithLength(object):
    """A generator with a length attribute"""
    def __init__(self, gen, length): self.gen, self.length = gen, length
    def __iter__(self): return self.gen
    def __len__(self): return self.length

################################################################################
def pretty_now():
    """Returns some thing like '2014-07-24 11:12:45 CEST+0200'"""
    now = datetime.datetime.now(dateutil.tz.tzlocal())
    return now.strftime("%Y-%m-%d %H:%M:%S %Z%z")