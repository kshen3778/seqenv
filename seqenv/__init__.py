b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.9.0'

# Built-in modules #
import multiprocessing
from collections import defaultdict

# Internal modules #
from seqenv.fasta import FASTA
from seqenv.seqsearch.parallel import ParallelSeqSearch
from seqenv.common.cache import property_cached

# Third party modules #
from Bio import Entrez
Entrez.email = "I don't know who will be running this script"

################################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.

    * `abundances`: If you have sample information, then you can provide the
                    'abundances' argument too. It should be a TSV file.

    * `N`: Use only the top `N` sequences in terms of their abundance, discard
           the rest. Only valid if abundances are provided.

    * 'seq_type': Either `nucl` or `prot`. Defaults to `nucl`.

    * `text_source`: Either `source` for isolation source terms or `abstract` for parsing
                     the publication abstracts.

    * `search_algo`: Either 'blast' or 'usearch' or ...

    * `search_db`: Either 'nt' or 'silva' or ...

    * `num_threads`: The number of threads. Default to the number of cores on the
                     current machine.

    * Sequence similarity search filtering options:
        - min_identity: Defaults to 0.97
        - e_value: Defaults to 0.0001
        - max_targets: Defaults to 10
        - min_coverage: Defaults to 0.97
    """

    def __init__(self, input_file,
                 abundances   = None,
                 N            = 1000,
                 seq_type     = 'nucl',
                 text_source  = 'source',
                 search_algo  = 'blast',
                 search_db    = 'nt',
                 num_threads  = None,
                 min_identity = 0.97,
                 e_value      = 0.0001,
                 max_targets  = 10,
                 min_coverage = 0.97,
                 ):
        # Base parameters #
        self.input_file = FASTA(input_file)
        self.abundances = abundances
        self.N = N
        self.seq_type = seq_type
        self.text_source = text_source
        # Search parameters #
        self.search_algo = search_algo
        self.search_db = search_db
        # Number of cores to use #
        if num_threads is None: self.num_threads = multiprocessing.cpu_count()/2
        else: self.num_threads = num_threads
        self.num_threads = min(self.num_threads, self.input_file.count)
        # Hit filtering parameters #
        self.min_identity = min_identity
        self.e_value      = e_value
        self.max_targets  = max_targets
        self.min_coverage = min_coverage

    @property_cached
    def orig_names_to_renamed(self):
        """A dictionary linking every sequence's name in the input FASTA to a
        new name following the scheme "C1", "C2", "C3" etc."""
        return {seq.id:"C%i"%i for i, seq in enumerate(self.input_file)}

    @property
    def renamed_fasta(self):
        """Make a new fasta file where every name in the input FASTA file is replaced
        with "C1", "C2", "C3" etc. Returns this new FASTA file."""
        renamed_fasta = FASTA(self.input_file.prefix_path + '_renamed.fasta')
        if renamed_fasta.exists: return renamed_fasta
        print "STEP 1: Generate mappings for sequence headers in FASTA file."
        self.input_file.rename_sequences(renamed_fasta, self.orig_names_to_renamed)
        return renamed_fasta

    @property
    def only_top_sequences(self):
        """Make a new fasta file where only the top N sequences are included
        (in terms of their abundance)."""
        if not self.abundances: return self.renamed_fasta
        only_top_fasta = FASTA(self.input_file.prefix_path + '_top.fasta')
        if only_top_fasta.exists: return only_top_fasta
        print "Using: " + self.renamed_fasta
        print "STEP 1B: Get the top %i sequences (in terms of their abundances)." % self.N
        #TODO
        return only_top_fasta

    @property_cached
    def filtering(self):
        """Return a dictionary with the filtering options for the sequence similarity
        search."""
        return {
            'min_identity': self.min_identity,
            'e_value':      self.e_value,
            'max_targets':  self.max_targets,
            'min_coverage': self.min_coverage,
        }

    @property_cached
    def search(self):
        """The similarity search object with all the relevant parameters."""
        return ParallelSeqSearch(input_fasta = self.only_top_sequences,
                                 seq_type    = self.seq_type,
                                 algorithm   = self.search_algo,
                                 database    = self.search_db,
                                 num_threads = self.num_threads,
                                 filtering   = self.filtering)

    @property
    def search_results(self):
        """For every sequence, search against a database and return the best hits
        after filtering."""
        # Check that the search was run #
        if not self.search.out_path.exists:
            print "Using: " + self.only_top_sequences
            print "STEP 2: Similarity search against the '%s' database" % self.search_db
            self.search.run()
            print "STEP 3: Filter out bad hits from the search results"
            self.search.filter()
        # Parse the results #
        return self.search.results

    @property_cached
    def seq_to_gis(self):
        """A dictionary linking every input sequence to a list of gi identifiers found."""
        result = defaultdict(list)
        for hit in self.search_results:
            seq_name = hit[0]
            gi = hit[1].split('|')[1]
            result[seq_name].append(gi)
        return result

    @property_cached
    def gi_to_text(self):
        """A dictionary linking every gi identifier to some unspecified text blob.
        Typically the text is its isolation source or a list of abstracts."""
        unique_gis = set(gi for gis in self.seq_to_gis.values() for gi in gis)
        result = {}
        for gi in unique_gis:
            gb_entry = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="xml")
            gb_record = Entrez.read(gb_entry)[0]
            gb_qualifiers = gb_record['GBSeq_feature-table'][0]['GBFeature_quals']
            for qualifier in gb_qualifiers:
                key, value = qualifier.items()
                1/0
                if key == 'isolation_source': result[gi] = value
        # TODO #
        if self.text_source == 'source': pass
        if self.text_source == 'abstract': pass

    @property_cached
    def gi_to_matches(self):
        """A dictionary linking every gi identifier to a list of
        regions of interest in the text (i.e. a match) and their meaning in terms of concepts.
        A 'concept' or 'entity' here would be an envo term such as 'ENVO:01000047'"""
        # Call the tagger
        t = tagger.Tagger()
        # Load the dictionary #
        t.LoadNames('data/envo_entities.tsv', 'data/envo_names.tsv')
        # Load a global blacklist #
        t.LoadGlobal('data/envo_global.tsv')
        # Tag all the text #
        # -27 is envo terms, -26 could be tissues, etc.
        # The second argument can be left empty in our case (per document blacklisting)
        l = t.GetMatches(python_string, "", [-27])
        start_pos, end_pos, concepts = l[0]
        type, id = concept[0] # -27, ENVO:01000047

    @property
    def gi_to_concepts(self, ):
        """Parse the matches data structure to extract the concepts and their respective counts"""
        # Loop over the counts
        pass
        # If backtracking is activated, add all the parent terms for every child term

    #--------------------------------------------------------------------------------------#
    @property
    def ncbi_results(self):
        """Using the search results, for every hit, download the relevant information from
        NCBI (e.g. abstract text)."""
        # Check that it was run #
        if not somefile.exists:
            print "Using %i results from the similarity search" % len(search_results)
            unique_ids = set(hit[1] for hit in self.search_results)
            # See envo_get_isolation_source.py
            # See envo_get_abstract.py
            print "STEP 4: Download data from NCBI."
        # Parse the results #
        return somestuff

    @property
    def tagger_results(self):
        """Using the NCBI results, for every piece of text, run the tagger on it."""
        # Check that it was run #
        if not somefile.exists:
            print "Using %i results from the NCBI data" % len(ncbi_results)
            print "STEP 5: Run the text mining tagger on the NCBI data."
        # Parse the results #
        return somestuff

    def generate_freq_tables(self):
        """Generate the frequencies tables..."""
        if not somefile.exists:
            print "Using %i results from the text mining results" % len(tagger_results)
            print "STEP 6: Generating main output."
        else:
            print "Final results already exist !"

    def generate_dot_files(self):
        """Generate the dot files for visualization in GraphViz..."""
        pass