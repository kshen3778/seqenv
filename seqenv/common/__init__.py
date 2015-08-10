# Built-in modules #
import tempfile

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