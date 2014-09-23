# Built-in modules #
import tempfile

################################################################################
def new_temp_path(**kwargs):
    handle = tempfile.NamedTemporaryFile(**kwargs)
    path = handle.name
    handle.close()
    return path