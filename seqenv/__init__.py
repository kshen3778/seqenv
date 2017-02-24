b'This module needs Python 2.7.x'

# Futures #
from __future__ import division

# Special variables #
__version__ = '1.3.0'
version_string = "seqenv version %s" % __version__

# Modules #
import sys, os
from seqenv.common.git import GitRepo

# Find the data dir #
module     = sys.modules[__name__]
module_dir = os.path.dirname(module.__file__) + '/'
repos_dir  = os.path.abspath(module_dir + '../') + '/'

# If we are in dev mode it's a git repo #
if os.path.exists(repos_dir + '.git/'): git_repo = GitRepo(repos_dir)
else:                                   git_repo = None

# Disable the X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Expose the main object #
from seqenv.analysis import Analysis