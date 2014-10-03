# Built-in modules #
import tempfile

# One liners #
flatter = lambda x: [item for sublist in x for item in sublist]

################################################################################
def new_temp_path(**kwargs):
    handle = tempfile.NamedTemporaryFile(**kwargs)
    path = handle.name
    handle.close()
    return path