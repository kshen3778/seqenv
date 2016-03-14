# Built-in modules #
import os, sh, glob, shutil

# Internal modules #
from seqenv.common.cache import property_cached

################################################################################
class DirectoryPath(str):

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.path)

    @classmethod
    def clean_path(cls, path):
        """Given a path, return a cleaned up version for initialization."""
        # Conserve 'None' object style #
        if path is None: return None
        # Don't nest DirectoryPaths or the like #
        if hasattr(path, 'path'): path = path.path
        # Expand the tilda #
        if "~" in path: path = os.path.expanduser(path)
        # Our standard is to end with a slash for directories #
        if not path.endswith('/'): path += '/'
        # Return the result #
        return path

    def __new__(cls, path, *args, **kwargs):
        """A DirectoryPath is in fact a string"""
        return str.__new__(cls, cls.clean_path(path))

    def __init__(self, path):
        self.path = self.clean_path(path)

    @property
    def name(self):
        """Just the directory name"""
        return os.path.basename(os.path.dirname(self.path))

    @property
    def exists(self):
        """Does it exist in the file system?"""
        return os.path.lexists(self.path) # Include broken symlinks

    def remove(self):
        if not self.exists: return False
        if self.is_symlink: return self.remove_when_symlink()
        shutil.rmtree(self.path, ignore_errors=True)
        return True

    def remove_when_symlink(self):
        if not self.exists: return False
        os.remove(self.path.rstrip('/'))
        return True

    def create(self, safe=False, inherit=True):
        # Create it #
        if not safe:
            os.makedirs(self.path)
            if inherit: os.chmod(self.path, self.directory.permissions.number)
        if safe:
            try:
                os.makedirs(self.path)
                if inherit: os.chmod(self.path, self.directory.permissions.number)
            except OSError: pass

################################################################################
class FilePath(str):
    """I can never remember all those darn `os.path` commands, so I made a class
    that wraps them with an easier and more pythonic syntax.

        path = FilePath('/home/root/text.txt')
        print path.extension
        print path.directory
        print path.filename

    You can find lots of the common things you would need to do with file paths.
    Such as: path.make_executable() etc etc."""

    def __new__(cls, path, *args, **kwargs):
        """A FilePath is in fact a string"""
        return str.__new__(cls, cls.clean_path(path))

    def __init__(self, path):
        self.path = self.clean_path(path)

    def __iter__(self): return open(self.path)

    def __len__(self):
        if self.path is None: return 0
        return self.count_lines

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
        """The directory containing this file"""
        if os.path.dirname(self.path) == "": return DirectoryPath(self)
        return DirectoryPath(os.path.dirname(self.path) + '/')

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

    def link_from(self, path, safe=False):
        """Make a link here pointing to another file somewhere else.
        The destination is hence self.path and the source is *path*."""
        if not safe:
            self.remove()
            return os.symlink(path, self.path)
        if safe:
            try: os.remove(self.path)
            except OSError: pass
            try: os.symlink(path, self.path)
            except OSError: pass

    def link_to(self, path, safe=False):
        """Create a link somewhere else pointing to this file.
        The destination is hence *path* and the source is self.path."""
        if not safe:
            if os.path.exists(path): os.remove(path)
            os.symlink(self.path, path)
        if safe:
            try: os.remove(path)
            except OSError: pass
            try: os.symlink(self.path, path)
            except OSError: pass
