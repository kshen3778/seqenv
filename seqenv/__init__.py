b'This module needs Python 2.7.x'

# Futures #
from __future__ import division

# Special variables #
__version__ = '1.0.4'
version_string = "seqenv version %s" % __version__

# Modules #
import sys, os

# Find the data dir #
module     = sys.modules[__name__]
module_dir = os.path.dirname(module.__file__) + '/'
repos_dir  = os.path.abspath(module_dir + '../') + '/'

# Expose the main object #
from seqenv.analysis import Analysis