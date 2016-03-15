# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing, warnings, unicodedata
from collections import defaultdict
import cPickle as pickle

# Internal modules #
from seqenv                    import module_dir, version_string, git_repo
from seqenv.fasta              import FASTA
from seqenv.outputs            import OutputGenerator
from seqenv.seqsearch.parallel import ParallelSeqSearch
from seqenv.common.cache       import property_cached
from seqenv.common.timer       import Timer
from seqenv.common.autopaths   import FilePath
from seqenv.common.database    import Database

# Third party modules #
import pandas

# Compiled modules #
try:
    import tagger
except ImportError:
    tagger = None
    msg = "Failed to import the tagger submodule, maybe it wasn't compiled correctly ?"
    warnings.warn(msg + " Seqenv won't work until that is fixed.", UserWarning)

################################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.
    The text below is somewhat redundant with what is in the README file and also
    the text in the command line script.

    * 'seq_type': Either `nucl` or `prot`. Defaults to `nucl`.

    * `search_algo`: Either 'blast' or 'usearch'. Defaults to `blast`.

    * `search_db`: The path to the database to search against. Defaults to `nt`.

    * `normalization`: Can be either of `flat`, `ui` or `upui`.
                        - If you choose `flat`, we will count every isolation source once,
                          even if the same text entry appears several time for the same inputs
                          sequence.
                        - If you choose `ui`, standing for unique isolation, we will uniquify
                          every frequency count depending on the text entry of its isolation
                          source.
                        - If you choose `upui`, standing for unique isolation and unique pubmed-ID,
                          we will uniquify the frequency counts based on the text entry of its
                          isolation source and the pubmed-ID from which the isolation text was
                          obtained.
                       This option defaults to `ui`.

    * `proportional`: Should we divide the counts of every input sequence by the
                      number of text entries that were associated to it.
                      Defaults to `True`.

    * `backtracking`: For every term identified by the tagger, we will propagate
                      frequency counts up the acyclic directed graph described by
                      the ontology. Defaults to `False`.

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

    def __repr__(self): return '<Analysis object on "%s" with %i sequences>' % \
                        (self.input_file.filename, self.input_file.count)

    def __init__(self, input_file,
                 seq_type      = 'nucl',
                 search_algo   = 'blast',
                 search_db     = 'nt',
                 normalization = 'ui',
                 proportional  = True,
                 backtracking  = False,
                 num_threads   = None,
                 out_dir       = None,
                 min_identity  = 0.97,
                 e_value       = 0.0001,
                 max_targets   = 10,
                 min_coverage  = 0.97,
                 abundances    = None,
                 N             = None):
        # Base parameters #
        self.input_file = FASTA(input_file)
        self.input_file.must_exist()
        # Abundance file #
        self.abundances = FilePath(abundances)
        if self.abundances: self.abundances.must_exist()
        # Other parameters #
        self.N = N
        self.seq_type = seq_type
        self.backtracking = bool(backtracking)
        self.proportional = bool(proportional)
        # Normalization parameters #
        options = ('flat', 'ui', 'upui')
        message = 'Normalization has to be one of %s' % (','.join(options))
        if normalization not in options: raise Exception(message)
        self.normalization = normalization
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
        # The object that can make the outputs for the user #
        self.outputs = OutputGenerator(self)

    def run(self):
        """A method to run the whole pipeline. As everything is coded in a functional
        style, we just need to make a call to `outputs.make_all` and everything will be
        generated automatically, in a reverse fashion."""
        print version_string + " (pid %i)" % os.getpid()
        print "The exact version of the code is: " + git_repo.short_hash
        self.timer.print_start()
        self.outputs.make_all()
        print "------------\nSuccess. Outputs are in '%s'" % self.out_dir
        self.timer.print_end()
        self.timer.print_total_elapsed()

    # --------------------------- In this section --------------------------- #
    # orig_names_to_renamed
    # renamed_to_orig
    # renamed_fasta
    # df_abundances
    # only_top_sequences
    # filtering
    # search
    # search_results

    @property_cached
    def orig_names_to_renamed(self):
        """A dictionary linking every sequence's name in the input FASTA to a
        new name following the scheme "C1", "C2", "C3" etc."""
        return {seq.id:"C%i"%i for i, seq in enumerate(self.input_file)}

    @property_cached
    def renamed_to_orig(self):
        """The opposite of the above dictionary"""
        return dict((v,k) for k,v in self.orig_names_to_renamed.items())

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

    @property_cached
    def df_abundances(self):
        """A pandas DataFrame object containing the abundance counts
        with OTUs as rows and sample names as columns."""
        assert self.abundances
        return pandas.io.parsers.read_csv(self.abundances, sep='\t', index_col=0, encoding='utf-8')

    @property
    def only_top_sequences(self):
        """Make a new fasta file where only the top N sequences are included
        (in terms of their abundance). Skipped if no abundance info is given."""
        if not self.abundances: return self.renamed_fasta
        only_top_fasta = FASTA(self.out_dir + 'top_seqs.fasta')
        if only_top_fasta.exists: return only_top_fasta
        # Use the default in case it was not provided #
        if self.N is None: N = 1000
        else:              N = int(self.N)
        # Print status #
        print "Using: " + self.renamed_fasta
        print "--> STEP 1B: Get the top %i sequences (in terms of their abundances)." % N
        # Check the user inputted value #
        if self.N is not None and N > self.input_file.count:
            msg = "You asked for the top %i sequences"
            msg += ", but your input file only contains %i sequences!"
            msg = msg % (self.N, self.input_file.count)
            warnings.warn(msg, UserWarning)
            N = self.input_file.count
        # Do it #
        ids = self.df_abundances.sum(axis=1).sort_values(ascending=False).index[0:N]
        ids = set([self.orig_names_to_renamed[x] for x in ids])
        self.renamed_fasta.extract_sequences(only_top_fasta, ids)
        self.timer.print_elapsed()
        return only_top_fasta

    @property_cached
    def filtering(self):
        """Return a dictionary with the filtering options for the sequence similarity
        search."""
        return {'min_identity': self.min_identity,
                'e_value':      self.e_value,
                'max_targets':  self.max_targets,
                'min_coverage': self.min_coverage}

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
            message = "--> STEP 2: Similarity search against the '%s' database with %i processes"
            print message % (self.search_db, self.num_threads)
            self.search.run()
            self.timer.print_elapsed()
            print "--> STEP 3: Filter out bad hits from the search results"
            self.search.filter()
            self.timer.print_elapsed()
            if self.search.out_path.count_bytes == 0:
                raise Exception("Found exactly zero hits after the similarity search.")
        # Parse the results #
        return self.search.results

    # --------------------------- In this section --------------------------- #
    # seq_to_gis
    # unique_gis
    # source_database
    # gis_with_text
    # unique_texts
    # text_to_matches
    # text_to_counts
    # seq_to_counts

    @property_cached
    def seq_to_gis(self):
        """A dictionary linking every input sequence to a list of gi identifiers found
        that are relating to it. If a sequence had no hits it links to an empty list.
        NB: You will get a KeyError if there are sequences in the search result files
        that are not present in the inputed fasta."""
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
        with open(seq_to_gis) as handle: return pickle.load(handle)

    @property_cached
    def unique_gis(self):
        """A set containing every GI that was found, considering all sequences combined."""
        return set(gi for gis in self.seq_to_gis.values() for gi in gis)

    @property_cached
    def source_db(self):
        """The sqlite3 database containing every GI number in all NCBI that has an
        isolation source associated to it. In addition, the pubmed-ID is listed too
        if there is one. If we don't have it locally already, we will go download it.
        Thus, the database containing two tables:

        CREATE TABLE "isolation"
            (
              "id"     INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
              "source" TEXT NOT NULL,
              "envos"  BLOB NOT NULL
            );
        CREATE TABLE "data"
            (
              "id"      INTEGER PRIMARY KEY NOT NULL,
              "isokey"  INTEGER NOT NULL REFERENCES "isolation"("id"),
              "pubmed"  INTEGER
            );
        CREATE INDEX "isolation_index" on "isolation" (source);
        CREATE INDEX 'main_index' on 'data' (id);"""
        path     = module_dir + 'data_sources/gi_db.sqlite3'
        drop_box = "ts5at7sISLFe9HxAyVjyemywNL78dMecTrNdoYmuD7DqSLUFxfpixCaPtvMZAOLB"
        retrieve = "https://dl.dropboxusercontent.com/content_link/%s/file?dl=1" % drop_box
        md5      = "e02d307bb099caced9f50843172fca26"
        database = Database(path, retrieve=retrieve, md5=md5, text_fact=bytes)
        self.timer.print_elapsed()
        return database

    @property_cached
    def gis_with_text(self):
        """A set containing every GI that was found and that had a source text associated."""
        return set(gi for gi in self.unique_gis if gi in self.source_db)

    @property_cached
    def unique_texts(self):
        """A set containing every unique isolation source text as unicode strings."""
        return set(self.source_db[gi][1] for gi in self.gis_with_text)

    @property_cached
    def text_to_matches(self):
        """A dictionary linking every isolation source found in this analysis
        to zero, one or several matches. You will not get a KeyError if you attempt to
        get the matches for an isolation source that had none."""
        text_to_matches = FilePath(self.out_dir + 'text_to_matches.pickle')
        # Check that it was run #
        if not text_to_matches.exists:
            print "--> STEP 5: Loading database with all NCBI isolation sources"
            print "Got %i unique GIs from search results" % len(self.unique_gis)
            print "Got %i unique GIs with an isolation source" % len(self.gis_with_text)
            print "Got %i unique isolation source texts" % len(self.unique_texts)
            self.timer.print_elapsed()
            print "--> STEP 6: Run the text mining tagger on all blobs."
            # Create the tagger #
            t = tagger.Tagger()
            # Load the dictionary #
            t.LoadNames(module_dir + 'data_envo/envo_entities.tsv',
                        module_dir + 'data_envo/envo_names.tsv')
            # Load a global blacklist #
            t.LoadGlobal(module_dir + 'data_envo/envo_global.tsv')
            # Tag all the texts, but the tagger only supports ascii #
            result = {}
            for text in self.unique_texts:
                ascii = unicodedata.normalize('NFKD', text).encode('ascii','ignore')
                result[text] = t.GetMatches(ascii, "", [-27])
            with open(text_to_matches, 'w') as handle: pickle.dump(result, handle)
            self.timer.print_elapsed()
            return result
        # Parse the results #
        with open(text_to_matches) as handle: return pickle.load(handle)

    @property_cached
    def text_to_counts(self):
        """A dictionary linking every isolation source to the concept counts.
        (dictionaries of concept:int). This is done by finding regions of interest in the text
        (i.e. a match). We can assign meaning to the matches in terms of concepts.
        A 'concept' or 'entity' here would be an envo term such as 'ENVO:01000047'
        You will get a KeyError if you attempt to get the counts for an isolation source
        that had none."""
        text_to_counts = FilePath(self.out_dir + 'text_to_counts.pickle')
        # Check that it was run #
        if not text_to_counts.exists:
            print "Got %i environmental term matches" % sum(map(len,self.text_to_matches.values()))
            print "--> STEP 7: Parsing the tagger results and counting terms."
            result = {}
            for text, matches in self.text_to_matches.items():
                if not matches: continue
                counts = defaultdict(float)
                for start_pos, end_pos, concepts in matches:
                    # This is the first place where the normalization technique used
                    # has to be thought about. For instance, at some point we had decided that
                    # every gi should adds up to 1.0 unless we have turned on backtracking.
                    # But this is not the case anymore in the current version!
                    ids = [concept_id for concept_type, concept_id in concepts]
                    score = 1 / len(ids) # Most of the time score is thus equal to 1
                    if self.backtracking: ids.extend([p for c in ids for p in self.child_to_parents[c]])
                    for concept_id in ids: counts[concept_id] += score
                result[text] = dict(counts)
            with open(text_to_counts, 'w') as handle: pickle.dump(result, handle)
            self.timer.print_elapsed()
            return result
        # Parse the results #
        with open(text_to_counts) as handle: return pickle.load(handle)

    @property_cached
    def seq_to_counts(self):
        """A dictionary linking every input sequence to its summed normalized concept
        counts dict, provided the input sequence had some hits, and at least one hit had
        a match. Otherwise it is empty.
        NB: What we want to account for is the fact that two GIs originating from the
        same sequence could be pointing to the same isolation source.
        In such case, we shall count the concepts from that isolation source only once.
        We can also avoid counting two GIs that are coming from the same study, if a pubmed
        number is available."""
        result = {}
        # Flat or Unique source are quite similar #
        if self.normalization == 'flat': set_or_list = list
        if self.normalization == 'ui':   set_or_list = set
        # It's just about using a list or a set in the right place #
        if self.normalization == 'flat' or self.normalization == 'ui':
            for seq, gis in self.seq_to_gis.items():
                texts = set_or_list(self.source_db[gi][1] for gi in gis if gi in self.source_db)
                counts = defaultdict(float)
                for text in texts:
                    if text not in self.text_to_counts: continue
                    for c,i in self.text_to_counts[text].items(): counts[c] += i
        # Unique source and unique pubmed #
        if self.normalization == 'uiup':
            for seq, gis in self.seq_to_gis.items():
                gis_w_text = [gi for gi in gis if gi in self.source_db]
                gi_to_tnp  = {gi: (self.source_db[gi][1], self.source_db[gi][2]) for gi in gis_w_text}
                texts      = set(text for text, pubmed in gi_to_tnp.values())
                find_first = lambda d,t: [(gi, tnp) for gi, tnp in d.values() if tnp[0]==t][0]
                gi_tnp_ut  = dict(find_first(gi_to_tnp, text) for text in texts)
                pubmeds    = set(pubmed for text, pubmed in gi_tnp_ut.values())
                find_first = lambda d,p: [(gi, tnp) for gi, tnp in d.values() if tnp[1]==p][0]
                gi_tnp_up  = dict(find_first(gi_to_tnp, pubmed) for pubmed in pubmeds)
                counts = defaultdict(float)
                for text, pubmed in gi_tnp_up:
                    if text not in self.text_to_counts: continue
                    for c,i in self.text_to_counts[text].items(): counts[c] += i
        # Proportional option #
        if self.proportional:
            tot_matches = sum([len(self.text_to_matches[text]) for text in texts])
            for k in counts: counts[k] /= tot_matches
        result[seq] = counts
        # Return #
        return result

    # --------------------------- In this section --------------------------- #
    # serial_to_concept
    # child_to_parents
    # concept_to_name

    @property_cached
    def serial_to_concept(self):
        """A dictionary linking every concept serial to its concept id.
        Every line in the file contains three columns: serial, concept_type, concept_id."""
        return dict(line.split()[0::2] for line in open(module_dir + 'data_envo/envo_entities.tsv'))

    @property_cached
    def child_to_parents(self):
        """A dictionary linking every concept id to a list of parent concept ids.
        Every line in the file contains two columns: child_serial, parent_serial"""
        result = defaultdict(list)
        with open(module_dir + 'data_envo/envo_groups.tsv') as handle:
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
        return dict(line.strip('\n').split('\t') for line in open(module_dir + 'data_envo/envo_preferred.tsv'))