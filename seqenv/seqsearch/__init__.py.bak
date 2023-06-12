# Built-in modules #
import multiprocessing

# Internal modules #
from seqenv.seqsearch.blast import BLASTquery
from seqenv.seqsearch.vsearch import VSEARCHquery
from seqenv.common.cache import property_cached
from seqenv.common.autopaths import FilePath

# Third party modules #

################################################################################
class SeqSearch(object):
    """A sequence similarity search. Could use different algorithms such
    as BLAST, VSEARCH, BLAT etc.

    Input: - List of sequences in a FASTA file
           - The type of the sequences
           - A database to search against
           - The type of algorithm to use
           - Number of threads to use
           - The desired output path
           - The filtering options:
             * BLAST supported:   - Minimum identity
                                  - E value
                                  - Maximum targets
                                  - Minimum query coverage (via manual output format)
             * VSEARCH supported: - ?
    Output: - Sorted list of identifiers in the database (object with significance value and identity attached)
    """

    def __init__(self, input_fasta, seq_type, database, algorithm='blast', num_threads=None, filtering=None, out_path=None):
        # Base parameters #
        self.input_fasta = input_fasta
        self.seq_type = seq_type
        self.database = database
        # Optional #
        self.algorithm = algorithm
        # The filtering options #
        if filtering is None: self.filtering = {}
        else: self.filtering = filtering
        # Number of cores to use #
        if num_threads is None: self.num_threads = multiprocessing.cpu_count()
        else: self.num_threads = num_threads
        # Output path #
        if out_path is None: self.out_path = FilePath(self.input_fasta.prefix_path + '.' + algorithm + 'out')
        else: self.out_path = FilePath(out_path)

    @property
    def query(self):
        """The similarity search object with all the relevant parameters."""
        if self.algorithm == 'blast': return self.blast_query
        if self.algorithm == 'vsearch': return self.vsearch_query
        raise NotImplemented(self.algorithm)

    @property_cached
    def blast_params(self):
        """A dictionary of options to pass to the blast executable.
        The params should depend on the filtering options."""
        # Initialize #
        params = {}
        # Defaults #
        params['-num_threads'] = self.num_threads
        params['-dust'] = 'no'
        # Conditionals #
        if 'e_value'      in self.filtering: params['-evalue']          = self.filtering['e_value']
        if 'max_targets'  in self.filtering: params['-max_target_seqs'] = self.filtering['max_targets']
        # Depends on sequence type #
        if self.seq_type == 'nucl':
            if 'min_identity' in self.filtering: params['-perc_identity'] = self.filtering['min_identity'] * 100
        # Output format #
        params['-outfmt'] = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids"
        # Return #
        return params

    @property_cached
    def blast_query(self):
        """Make a BLAST search object."""
        # Sequence type #
        if self.seq_type == 'nucl': blast_algo = 'blastn'
        if self.seq_type == 'prot': blast_algo = 'blastp'
        # The query object #
        query = BLASTquery(self.input_fasta, self.database, self.blast_params, blast_algo, "plus", self.out_path)
        return query

    @property_cached
    def vsearch_query(self):
        """Make a VSEARCH search object."""
        params = {}
        query = VSEARCHquery(self.input_fasta, self.database, params)
        return query

    def run(self):
        """Run the search"""
        return self.query.run()

    def filter(self):
        """Filter the results accordingly"""
        return self.query.filter(self.filtering)

    @property
    def results(self):
        """Parse the results."""
        return self.query.results