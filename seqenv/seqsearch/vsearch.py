# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from __future__ import absolute_import
from seqenv.common.autopaths import FilePath
from seqenv.fasta import FASTA

# Third party modules #
import sh
from six.moves import map

###############################################################################
class VSEARCHquery(object):
    """A vsearch job."""

    def __init__(self, query_path, db_path,
                 params = None,
                 out_path = None,
                 executable = None):
        # Save attributes #
        self.query = FASTA(query_path)
        self.db = db_path
        self.params = params if params else {}
        self.executable = FilePath(executable)
        # Output #
        if out_path is None: self.out_path = self.query.prefix_path + '.vsearchout'
        elif out_path.endswith('/'): self.out_path = out_path + self.query.prefix + '.vsearchout'
        else: self.out_path = out_path
        self.out_path = FilePath(self.out_path)

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        # Return #
        return list(map(str, cmd))

    def run(self):
        sh.Command(self.command[0])(self.command[1:])