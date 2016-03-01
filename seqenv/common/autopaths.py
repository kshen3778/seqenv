# Built-in modules #
import os, sh, glob

# Internal modules #
from seqenv.common.cache import property_cached

################################################################################
class FilePath(str):
    """A file path somewhere in the file system. Useful methods for dealing
    with such paths are included. For instance, you can ask for `path.extension`"""

    @classmethod
    def clean_path(cls, path):
        """Given a path, return a cleaned up version for initialization"""
        # Conserve None object style #
        if path is None: return None
        # Don't nest FilePaths or the like #
        if hasattr(path, 'path'): path = path.path
        # Expand tilda #
        if "~" in path: path = os.path.expanduser(path)
        # Expand star #
        if "*" in path:
            matches = glob.glob(path)
            if len(matches) < 1: raise Exception("Found exactly no files matching '%s'" % path)
            if len(matches) > 1: raise Exception("Found several files matching '%s'" % path)
            path = matches[0]
        # Return the result #
        return path

    def __new__(cls, path, *args, **kwargs):
        """A FilePath is in fact a string"""
        return str.__new__(cls, cls.clean_path(path))

    def __init__(self, path):
        self.path = self.clean_path(path)

    def __iter__(self): return open(self.path)
    def __len__(self):
        if self.path is None: return 0
        return self.count_lines

    @property_cached
    def count_lines(self):
        return int(sh.wc('-l', self.path).split()[0])

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
        if os.path.dirname(self.path) == "": return './'
        return os.path.dirname(self.path) + '/'

    @property
    def count_bytes(self):
        """The number of bytes"""
        if not self.exists: return 0
        return os.path.getsize(self.path)

    def remove(self):
        if not self.exists: return False
        os.remove(self.path)
        return True

    def write(self, contents):
        with open(self.path, 'w') as handle: handle.write(contents)

    def writelines(self, contents):
        with open(self.path, 'w') as handle: handle.writelines(contents)

    def must_exist(self):
        """Raise an exception if the path doesn't exist."""
        if not self.exists: raise Exception("The file path '%s' does not exist." % self.path)