# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing

# Internal modules #
from seqenv.common.autopaths import FilePath
from seqenv.fasta import FASTA

# Third party modules #
import sh

###############################################################################
class BLASTquery(object):
    """A blast job. Possibly the standard BLAST algorithm or
       BLASTP or BLASTX etc. Typically you could use it like this:

            import sys, os
            records_path = os.path.expanduser(sys.argv[1])
            centers_path = 'centers.fasta'
            db = parallelblast.BLASTdb(centers_path)
            db.makeblastdb()
            params = {'executable': "~/share/blastplus/blastn",
                      '-outfmt': 0,
                      '-evalue': 1e-2,
                      '-perc_identity': 97,
                      '-num_threads': 16}
            search = parallelblast.BLASTquery(records_path, db, params)
            search.run()
       """

    def __repr__(self): return '<%s on "%s">' % (self.__class__.__name__, self.query_path)

    def __init__(self, query_path, db_path,
                 params = None,
                 algorithm = "blastn",
                 version = "plus" or "legacy",
                 out_path = None,
                 executable = None):
        # Save attributes #
        self.query = FASTA(query_path)
        self.db = FilePath(db_path)
        self.version = version
        self.algorithm = algorithm
        self.params = params if params else {}
        self.executable = FilePath(executable)
        # Output #
        if out_path is None: self.out_path = self.query.prefix_path + '.blastout'
        elif out_path.endswith('/'): self.out_path = out_path + self.query.prefix + '.blastout'
        else: self.out_path = out_path
        self.out_path = FilePath(self.out_path)
        # Defaults #
        cpus = multiprocessing.cpu_count()
        if self.version == 'plus':
            if '-num_threads' not in self.params: self.params['-num_threads'] = cpus
        if self.version == 'legacy':
            if '-a' not in self.params: self.params['-a'] = cpus

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        elif self.version == 'legacy': cmd = ["blastall", '-p', self.algorithm]
        else: cmd = [self.algorithm]
        # Essentials #
        if self.version == 'legacy': cmd += ['-d',  self.db, '-i',     self.query, '-o',   self.out_path]
        if self.version == 'plus':   cmd += ['-db', self.db, '-query', self.query, '-out', self.out_path]
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return map(str, cmd)

    def run(self):
        sh.Command(self.command[0])(self.command[1:])
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")

    def run_parallel(self):
        """Special method to use GNU parallel to run the query"""
        self.params['-num_threads'] = 1
        self.query = "{}"
        # --block $sizeChunksString --recstart '>' --pipe
        sh.parallel(self.query.cmd)