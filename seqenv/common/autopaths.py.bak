# Built-in modules #
import os, sh, glob, shutil, stat

# Internal modules #
from seqenv.common import unzip
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

    @property
    def absolute_path(self):
        """The absolute path starting with a `/`"""
        return os.path.abspath(self.path) + '/'

    @property
    def is_symlink(self):
        """Is this directory a symbolic link to an other directory?"""
        return os.path.islink(self.path.rstrip('/'))

    @property
    def directory(self):
        """The full path of directory containing this one"""
        return DirectoryPath(os.path.dirname(os.path.dirname(self.path)))

    @property
    def permissions(self):
        """Convenience object for dealing with permissions."""
        return FilePermissions(self.path)

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

    def create_if_not_exists(self):
        if not self.exists: self.create()

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
    def absolute_path(self):
        """The absolute path starting with a `/`"""
        return os.path.abspath(self.path)

    @property
    def directory(self):
        """The directory containing this file"""
        # The built-in function #
        directory = os.path.dirname(self.path)
        # Maybe we need to go the absolute path way #
        if not directory: directory = os.path.dirname(self.absolute_path)
        # Return #
        return DirectoryPath(directory+'/')

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

    @property
    def size(self):
        """Human readable file size"""
        return Filesize(self.count_bytes)

    def link_from(self, path, safe=False):
        """Make a link here pointing to another file somewhere else.
        The destination is hence self.path and the source is *path*."""
        # Get source and destination #
        source      = path
        destination = self.path
        # Do it #
        if not safe:
            if os.path.exists(destination): os.remove(destination)
            os.symlink(source, destination)
        # Do it safely #
        if safe:
            try: os.remove(destination)
            except OSError: pass
            try: os.symlink(source, destination)
            except OSError: pass

    def link_to(self, path, safe=False, absolute=True):
        """Create a link somewhere else pointing to this file.
        The destination is hence *path* and the source is self.path."""
        # Get source and destination #
        if absolute: source = self.absolute_path
        else:        source = self.path
        destination = path
        # Do it #
        if not safe:
            if os.path.exists(destination): os.remove(destination)
            os.symlink(source, destination)
        # Do it safely #
        if safe:
            try: os.remove(destination)
            except OSError: pass
            try: os.symlink(source, destination)
            except OSError: pass

    def unzip_to(self, destination=None, inplace=False):
       """Make an unzipped version of the file at a given path"""
       return unzip(self.path, destination=destination, inplace=inplace)

################################################################################
class FilePermissions(object):
    """Container for reading and setting a files permissions"""

    def __init__(self, path):
        self.path = path

    @property
    def number(self):
        """The permission bits as an octal integer"""
        return os.stat(self.path).st_mode & 0777

    def make_executable(self):
        return os.chmod(self.path, os.stat(self.path).st_mode | stat.S_IEXEC)

    def only_readable(self):
        """Remove all writing privileges"""
        return os.chmod(self.path, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)

################################################################################
class Filesize(object):
    """
    Container for a size in bytes with a human readable representation
    Use it like this:

        >>> size = Filesize(123123123)
        >>> print size
        '117.4 MiB'
    """

    chunk      = 1000 # Could be 1024 if you like old-style sizes
    units      = ['bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB']
    precisions = [0, 0, 1, 2, 2, 2]

    def __init__(self, size):
        self.size = size

    def __int__(self):
        return self.size

    def __eq__(self, other):
        return self.size == other

    def __str__(self):
        if self.size == 0: return '0 bytes'
        from math import log
        unit = self.units[min(int(log(self.size, self.chunk)), len(self.units) - 1)]
        return self.format(unit)

    def format(self, unit):
        # Input checking #
        if unit not in self.units: raise Exception("Not a valid file size unit: %s" % unit)
        # Special no plural case #
        if self.size == 1 and unit == 'bytes': return '1 byte'
        # Compute #
        exponent      = self.units.index(unit)
        quotient      = float(self.size) / self.chunk**exponent
        precision     = self.precisions[exponent]
        format_string = '{:.%sf} {}' % (precision)
        # Return a string #
        return format_string.format(quotient, unit)