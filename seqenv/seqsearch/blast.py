# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing, threading, shutil

# Internal modules #
from seqenv.fasta import FASTA
from seqenv.common import new_temp_path
from seqenv.common.autopaths import FilePath

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

       You can also call search.non_block_run() to run maybe searches in parallel.
       """

    def __init__(self, query_path, db_path,
                 params = None,
                 algorithm = "blastn",
                 version = "plus" or "legacy",
                 out_path = None,
                 executable = None):
        # Save attributes #
        self.path = query_path
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
        self.cpus = multiprocessing.cpu_count()
        if self.version == 'plus':
            if '-num_threads' not in self.params: self.params['-num_threads'] = self.cpus
        if self.version == 'legacy':
            if '-a' not in self.params: self.params['-a'] = self.cpus

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

    def non_block_run(self):
        """Special method to run the query in a thread without blocking."""
        self.thread = threading.Thread(target=self.run)
        self.thread.start()

    def wait(self):
        """If you have run the query in a non-blocking way, call this method to pause
        until the query is finished."""
        self.thread.join()

    def filter(self, filtering):
        """We can do some special filtering on the results. For the moment only minimum coverage."""
        # Conditions #
        if 'min_coverage' not in filtering: return
        if 'qcovs' not in self.params['-outfmt']: raise Exception()
        # Iterator #
        def filter_lines(blastout):
            threshold = filtering['min_coverage'] * 100
            position = self.params['-outfmt'].strip('"').split().index('qcovs') - 1
            for line in blastout:
                coverage = line.split()[position]
                if coverage < threshold: continue
                else: yield line
        # Do it #
        temp_path = new_temp_path()
        with open(temp_path, 'w') as handle: handle.writelines(filter_lines(self.out_path))
        os.remove(self.out_path)
        shutil.move(temp_path, self.out_path)

    @property
    def results(self):
        """Parse the results."""
        for line in self.out_path: yield line.split()
