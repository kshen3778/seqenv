# Built-in modules #
import multiprocessing

# Internal modules #
from seqenv.seqsearch.blast import BLASTquery
from seqenv.seqsearch.usearch import USEARCHquery
from seqenv.common.cache import property_cached

# Third party modules #

################################################################################
class SeqSearch(object):
    """A sequence similarity search. Could use different algorithms such
    as BLAST, USEARCH, BLAT etc.

    Operates by chopping in the input up into smaller pieces and running the
    algorithm on each piece separately, finally joining the outputs.

    Input: - List of sequences
           - The type of the sequences
           - A database to search against
           - The type of algorithm to use
           - Number of threads to use
           - The filtering options:
             * BLAST supported:   - Minimum identity
                                  - E value
                                  - Maximum targets
             * USEARCH supported: - ?
             * General: Minimum query coverage
    Output: - Sorted list of identifiers in the database (object with significance value and identity attached)

    Other ideas:
    @property
    def word_size(self):
        #Depends on the percent identity
        pass
    """

    def __init__(self, input_fasta, seq_type, database, algorithm='blast', num_threads=None, filtering=None):
        # Base parameters #
        self.input_fasta = input_fasta
        self.seq_type = seq_type
        self.database = database
        self.algorithm = algorithm
        self.filtering = filtering
        # Number of cores to use #
        if num_threads is None: self.num_threads = multiprocessing.cpu_count()
        else: self.num_threads = num_threads

    @property
    def query(self):
        """The similarity search object with all the relevant parameters."""
        if self.algorithm == 'blast': return self.blast_query
        if self.algorithm == 'usearch': return self.usearch_query
        raise NotImplemented(self.algorithm)

    @property_cached
    def blast_query(self):
        """Make a BLAST search object."""
        # The params should depend on the filtering options #
        params = {'-dust': 'no', '-outfmt': '6'}
        # Sequence type #
        if self.seq_type == 'nucl': blast_algo = 'blastn'
        if self.seq_type == 'prot': blast_algo = 'blastp'
        # The query object #
        query = BLASTquery(self.input_fasta, self.database, params, blast_algo, version="plus")
        return query

    @property_cached
    def usearch_query(self):
        """Make a USEARCH search object."""
        params = {}
        query = USEARCHquery(self.input_fasta, self.database, params)
        return query

    def run(self):
        """Run the search"""
        # Parallel version #
        if self.num_threads > 1: return self.run_parallel()
        # Single threaded version #
        return self.query.run()

    def run_parallel(self):
        """Use GNU parallel to run the search job"""
        return self.query.run_parallel()

    def filter(self):
        raise NotImplemented('')