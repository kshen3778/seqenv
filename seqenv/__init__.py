b'This module needs Python 2.7.x'

# Futures #
from __future__ import division

# Special variables #
__version__ = '1.0.1'
version_string = "seqenv version %s" % __version__

# Built-in modules #
import os, inspect, multiprocessing, gzip
from collections import defaultdict
import cPickle as pickle

# Internal modules #
from seqenv.fasta import FASTA
from seqenv.outputs import OutputGenerator
from seqenv.seqsearch.parallel import ParallelSeqSearch
from seqenv.common.cache import property_cached
from seqenv.common.timer import Timer
from seqenv.common.autopaths import FilePath

# Compiled modules #
import tagger

# Third party modules #
import pandas
from tqdm import tqdm

# Constants #
total_gis_with_source = 13658791

# Find the data dir #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
possible = current_dir + 'data/envo_preferred.tsv'
if os.path.exists(possible): data_dir = current_dir + 'data/'
else: possible = '' #TODO

################################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.
    The text below is somewhat redundant with what is in the README file and in
    the command line script.

    * 'seq_type': Either `nucl` or `prot`. Defaults to `nucl`.

    * `search_algo`: Either 'blast' or 'usearch'. Defaults to `blast`.

    * `search_db`: The path to the database to search against. Defaults to `nt`.

    * `backtracking`: For every term identified by the tagger, we will propagate
                      frequency counts up the acyclic directed graph described by
                      the ontology. Defaults to `False`.

    * `normalization`: Should we divide the counts of every input sequence by the
                       number of text entries that were associated to it. Defaults to `True`.

    * `num_threads`: The number of threads. Default to the number of cores on the
                     current machine.

    * `out_dir`: Place all the outputs in the specified directory.
                 Defaults to the input file's directory.

    * Sequence similarity search filtering options:
        - min_identity: Defaults to 0.97
        - e_value: Defaults to 0.0001
        - max_targets: Defaults to 10
        - min_coverage: Defaults to 0.97

    * `abundances`: If you have sample information, then you can provide the
                    'abundances' argument too. It should be a TSV file with OTUs
                    as rows and sample names as columns.

    * `N`: Use only the top `N` sequences in terms of their abundance, discard
           the rest. Only valid if abundances are provided.
    """

    def __init__(self, input_file,
                 seq_type      = 'nucl',
                 search_algo   = 'blast',
                 search_db     = 'nt',
                 backtracking  = False,
                 normalization = True,
                 num_threads   = None,
                 out_dir       = None,
                 min_identity  = 0.97,
                 e_value       = 0.0001,
                 max_targets   = 10,
                 min_coverage  = 0.97,
                 abundances    = None,
                 N             = 1000):
        # Base parameters #
        self.input_file = FASTA(input_file)
        self.input_file.must_exist()
        # Abundance file #
        self.abundances = FilePath(abundances)
        if self.abundances: self.abundances.must_exist()
        # Other parameters #
        self.N = int(N)
        self.seq_type = seq_type
        self.backtracking = bool(backtracking)
        self.normalization = bool(normalization)
        # Search parameters #
        self.search_algo = search_algo
        self.search_db = search_db
        # Number of cores to use #
        if num_threads is None: self.num_threads = min(multiprocessing.cpu_count(), 32)
        else: self.num_threads = int(num_threads)
        self.num_threads = min(self.num_threads, self.input_file.count)
        # Hit filtering parameters #
        self.min_identity = float(min_identity)
        self.e_value      = float(e_value)
        self.max_targets  = int(max_targets)
        self.min_coverage = float(min_coverage)
        # Time the pipeline execution #
        self.timer = Timer()
        # Keep all outputs in a directory #
        if out_dir is None: self.out_dir = self.input_file.directory
        else: self.out_dir = out_dir
        if not self.out_dir.endswith('/'): self.out_dir += '/'
        if not os.path.exists(self.out_dir): os.makedirs(self.out_dir)
        # The object that can make the outputs #
        self.outputs = OutputGenerator(self)

    def run(self):
        """A method to run the whole pipeline. As everything is coded in a functional
        style, we just need to make a call to `outputs.make_all` and everything will
        generated automatically"""
        print version_string + " (pid %i)" % os.getpid()
        self.timer.print_start()
        self.outputs.make_all()
        print "------------\nSuccess. Outputs are in '%s'" % self.out_dir
        self.timer.print_end()
        self.timer.print_total_elapsed()

    #-------------------------------------------------------------------------#
    @property_cached
    def orig_names_to_renamed(self):
        """A dictionary linking every sequence's name in the input FASTA to a
        new name following the scheme "C1", "C2", "C3" etc."""
        return {seq.id:"C%i"%i for i, seq in enumerate(self.input_file)}
    @property_cached
    def renamed_to_orig(self): return dict((v,k) for k,v in self.orig_names_to_renamed.items())

    @property
    def renamed_fasta(self):
        """Make a new fasta file where every name in the input FASTA file is replaced
        with "C1", "C2", "C3" etc. Returns this new FASTA file."""
        renamed_fasta = FASTA(self.out_dir + 'renamed.fasta')
        if renamed_fasta.exists: return renamed_fasta
        print "--> STEP 1: Parse the input FASTA file."
        self.input_file.rename_sequences(renamed_fasta, self.orig_names_to_renamed)
        self.timer.print_elapsed()
        return renamed_fasta

    @property
    def only_top_sequences(self):
        """Make a new fasta file where only the top N sequences are included
        (in terms of their abundance)."""
        if not self.abundances: return self.renamed_fasta
        only_top_fasta = FASTA(self.out_dir + 'top_seqs.fasta')
        if only_top_fasta.exists: return only_top_fasta
        print "Using: " + self.renamed_fasta
        print "--> STEP 1B: Get the top %i sequences (in terms of their abundances)." % self.N
        ids = self.df_abundances.sum(axis=1).sort(inplace=False, ascending=False).index[0:self.N]
        ids = set([self.orig_names_to_renamed[x] for x in ids])
        self.renamed_fasta.extract_sequences(only_top_fasta, ids)
        self.timer.print_elapsed()
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
            print "--> STEP 2: Similarity search against the '%s' database with %i processes" % (self.search_db, self.num_threads)
            self.search.run()
            self.timer.print_elapsed()
            print "--> STEP 3: Filter out bad hits from the search results"
            self.search.filter()
            self.timer.print_elapsed()
            if self.search.out_path.count_bytes == 0:
                raise Exception("Found exactly zero hits after the similarity search.")
        # Parse the results #
        return self.search.results

    #-------------------------------------------------------------------------#
    @property_cached
    def seq_to_gis(self):
        """A dictionary linking every input sequence to a list of gi identifiers found
        that are relating to it. You will get a KeyError if there are sequences in the search
        result files that are not present in the inputed fasta."""
        seq_to_gis = FilePath(self.out_dir + 'seq_to_gis.pickle')
        # Check that is was run #
        if not seq_to_gis.exists:
            self.search_results
            print "--> STEP 4: Parsing the search results"
            result = {seq:[] for seq in self.only_top_sequences.ids}
            for hit in self.search_results:
                seq_name = hit[0]
                gi = hit[1].split('|')[1]
                result[seq_name].append(gi)
            with open(seq_to_gis, 'w') as handle: pickle.dump(result, handle)
            self.timer.print_elapsed()
            return result
        # Parse the results #
        with open(seq_to_gis, 'r') as handle: return pickle.load(handle)

    @property_cached
    def gi_to_text(self):
        """A dictionary linking every gi identifier in NCBI to its isolation sourcet
        test, provided it has one."""
        print "--> STEP 5: Loading all NCBI isolation sources in RAM"
        result = {}
        with gzip.open(data_dir + 'gi_to_source.tsv.gz') as handle:
            for i in tqdm(xrange(total_gis_with_source), total=total_gis_with_source):
                gi, source = handle.next().lstrip('GI:').rstrip('\n').split('\t')
                result[gi] = source
        self.timer.print_elapsed()
        return result

    @property_cached
    def gi_to_matches(self):
        """A dictionary linking every gi to zero, one or several matches.
        When you call `t.GetMatches(text, "", [-27])` you get a list back.
        The second argument can be left empty in our case (per document blacklisting)
        The result is something like:
        - [(53, 62, ((-27, 'ENVO:00002001'),)), (64, 69, ((-27, 'ENVO:00002044'),))]
        The number -27 is ENVO terms, -26 could be tissues, etc.
        """
        gi_to_matches = FilePath(self.out_dir + 'gi_to_matches.pickle')
        # Check that it was run #
        if not gi_to_matches.exists:
            unique_gis = set(gi for gis in self.seq_to_gis.values() for gi in gis)
            print "Got %i GIs from search results" % len(unique_gis)
            print "Got %i GIs with isolation source" % len(self.gi_to_text)
            print "--> STEP 6: Run the text mining tagger on all blobs."
            # Create the tagger #
            t = tagger.Tagger()
            # Load the dictionary #
            t.LoadNames(data_dir + 'envo_entities.tsv', data_dir + 'envo_names.tsv')
            # Load a global blacklist #
            t.LoadGlobal(data_dir + 'envo_global.tsv')
            # Tag all the text #
            result = {}
            for gi in unique_gis:
                text = self.gi_to_text.get(gi)
                if text is None:
                    result[gi] = []
                    continue
                result[gi] = t.GetMatches(text, "", [-27])
            with open(gi_to_matches, 'w') as handle: pickle.dump(result, handle)
            self.timer.print_elapsed()
            return result
        # Parse the results #
        with open(gi_to_matches, 'r') as handle: return pickle.load(handle)

    @property_cached
    def gi_to_counts(self):
        """A dictionary linking every `gi` identifier to the concept counts.
        (dictionaries of concept:int). This is done by finding regions of interest in the text
        (i.e. a match). We can assign meaning to the matches in terms of concepts.
        A 'concept' or 'entity' here would be an envo term such as 'ENVO:01000047'"""
        gi_to_counts = FilePath(self.out_dir + 'gi_to_counts.pickle')
        # Check that it was run #
        if not gi_to_counts.exists:
            print "Using matches from %i gi entries" % len(self.gi_to_matches)
            print "--> STEP 7: Parsing the tagger results and counting terms."
            result = {}
            for gi, matches in self.gi_to_matches.items():
                if not matches: continue
                counts = defaultdict(int)
                for start_pos, end_pos, concepts in matches:
                    ids = [concept_id for concept_type, concept_id in concepts]
                    score = 1 / len(ids) # Every gi adds up to one unless we have backtracking
                    if self.backtracking: ids.extend([p for c in ids for p in self.child_to_parents[c]])
                    for concept_id in ids: counts[concept_id] += score
                result[gi] = dict(counts)
            with open(gi_to_counts, 'w') as handle: pickle.dump(result, handle)
            self.timer.print_elapsed()
            return result
        # Parse the results #
        with open(gi_to_counts, 'r') as handle: return pickle.load(handle)

    @property_cached
    def seq_to_counts(self):
        """A dictionary linking every input sequence to its summed normalized concept counts dict,
        provided the input sequence had some hits, and a hit had a match, otherwise it is empty."""
        result = {}
        for seq, gis in self.seq_to_gis.items():
            counts = defaultdict(int)
            for gi in gis:
                if gi not in self.gi_to_counts: continue
                for c,i in self.gi_to_counts[gi].items(): counts[c] += i
            if self.normalization:
                tot_matches = sum([len(self.gi_to_matches[gi]) for gi in gis])
                for k in counts: counts[k] /= tot_matches
            result[seq] = counts
        return result

    #-------------------------------------------------------------------------#
    @property_cached
    def serial_to_concept(self):
        """A dictionary linking every concept serial to its concept id.
        Every line in the file contains three columns: serial, concept_type, concept_id
        This could possibly overflow the memory when we come with NCBI taxonomy etc."""
        return dict(line.split()[0::2] for line in open(data_dir + 'envo_entities.tsv'))

    @property_cached
    def child_to_parents(self):
        """A dictionary linking every concept id to a list of parent concept ids.
        Every line in the file contains two columns: child_serial, parent_serial"""
        result = defaultdict(list)
        with open(data_dir + 'envo_groups.tsv') as handle:
            for line in handle:
                child_serial, parent_serial = line.split()
                child_concept = self.serial_to_concept[child_serial]
                parent_concept = self.serial_to_concept[parent_serial]
                result[child_concept].append(parent_concept)
        return result

    @property_cached
    def concept_to_name(self):
        """A dictionary linking the concept id to relevant names. In this case ENVO terms.
        Hence, ENVO:00000095 would be linked to 'lava field'"""
        return dict(line.strip('\n').split('\t') for line in open(data_dir + 'envo_preferred.tsv'))

    @property_cached
    def df_abundances(self):
        """A pandas DataFrame object containing the abundance counts
        with OTUs as rows and sample names as columns."""
        assert self.abundances
        return pandas.io.parsers.read_csv(self.abundances, sep='\t', index_col=0, encoding='utf-8')