# Built-in modules #
import os

################################################################################
class FilePath(str):
    """A file path somewhere in the file system. Useful methods for dealing
    with such paths are included. For instance, you can ask for `path.extension`"""

    def __new__(cls, path, *args, **kwargs):
        if path is None: return None
        if isinstance(path, FilePath): path = path.path
        return str.__new__(cls, path)

    def __init__(self, path):
        if isinstance(path, FilePath): path = path.path
        self.path = path

    def __iter__(self): return open(self.path)

    @property
    def exists(self):
        """Does it exist in the file system. Returns True or False."""
        return os.path.lexists(self.path)

    @property
    def prefix_path(self):
        """The full path without the (last) extension and trailing period"""
        return str(os.path.splitext(self.path)[0])

    @property
    def prefix(self):
        """Just the filename without the (last) extension and trailing period"""
        return str(os.path.basename(self.prefix_path))

    @property
    def filename(self):
        """Just the filename with the extension"""
        return str(os.path.basename(self.path))

    @property
    def extension(self):
        """The extension with the leading period"""
        return os.path.splitext(self.path)[1]

    @property
    def directory(self):
        """The containing directory"""
        return os.path.dirname(self.path) + '/'

    @property
    def count_bytes(self):
        """The number of bytes"""
        if not self.exists: return 0
        return os.path.getsize(self.path)

