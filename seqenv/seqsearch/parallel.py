# Built-in modules #
import subprocess

# Internal modules #
from seqenv.seqsearch import SeqSearch
from seqenv.seqsearch.blast import BLASTquery
from seqenv.seqsearch.vsearch import VSEARCHquery
from seqenv.common.cache import property_cached
from seqenv.fasta.splitable import SplitableFASTA

# Third party modules #

################################################################################
class ParallelSeqSearch(SeqSearch):
    """The same thing as a SeqSearch but operates by chopping the input up into
    smaller pieces and running the algorithm on each piece separately, finally
    joining the outputs.
    In addition, the pieces can be run separately on the local machine, or distributed to different
    compute nodes using the SLURM system.
    """

    @property_cached
    def splitable(self):
        """The input fasta file, but with the ability to split it."""
        return SplitableFASTA(self.input_fasta, self.num_threads)

    @property
    def queries(self):
        """A list of all the queries to run."""
        if self.algorithm == 'blast': return self.blast_queries
        if self.algorithm == 'usearch': return self.vsearch_queries
        raise NotImplemented(self.algorithm)

    @property_cached
    def blast_queries(self):
        """Make all BLAST search objects."""
        # Only use one CPU per query #
        params = self.blast_params.copy()
        params['-num_threads'] = 1
        # Sequence type #
        if self.seq_type == 'nucl': blast_algo = 'blastn'
        if self.seq_type == 'prot': blast_algo = 'blastp'
        # Make many queries #
        return [BLASTquery(p, self.database, params, blast_algo, version="plus") for p in self.splitable.parts]

    @property_cached
    def vsearch_queries(self):
        """Make all VSEARCH search objects."""
        # The params should depend on the filtering options #
        params = {}
        # Make many queries #
        return [VSEARCHquery(p, self.database, params) for p in self.splitable.parts]

    def run(self):
        """Run the search"""
        self.splitable.split()
        for query in self.queries: query.non_block_run()
        for query in self.queries: query.wait()
        self.join_outputs()

    def join_outputs(self):
        """Join the outputs"""
        command = "cat %s > %s" % (' '.join(q.out_path for q in self.queries), self.out_path)
        subprocess.call(command, shell=True)