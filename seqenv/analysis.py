# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing, warnings, marshal
from collections import defaultdict
import cPickle as pickle

# Internal modules #
from seqenv                    import module_dir, version_string, git_repo
from seqenv.fasta              import FASTA
from seqenv.ontology           import Ontology
from seqenv.outputs            import OutputGenerator
from seqenv.seqsearch.parallel import ParallelSeqSearch
from seqenv.common.cache       import property_cached, cached, class_property
from seqenv.common.timer       import Timer
from seqenv.common.autopaths   import FilePath
from seqenv.common.database    import Database

# Third party modules #
import pandas

################################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.
    The text below is somewhat redundant with what is in the README file and also
    the text in the command line script.

    * 'seq_type': Either `nucl` or `prot`. Defaults to `nucl`.

    * `search_algo`: Either 'blast' or 'usearch'. Defaults to `blast`.

    * `search_db`: The path to the database to search against. Defaults to `nt`.

    * `normalization`: Can be either of `flat`, `ui` or `upui`.
      This option defaults to `flat`.
         * If you choose `flat`, we will count every isolation source independently,
           even if the same text appears several times for the same input sequence.
         * If you choose `ui`, standing for unique isolation, we will count every
           identical isolation source only once within the same input sequence.
         * If you choose `upui`, standing for unique isolation and unique pubmed-ID,
           we will uniquify the counts based on the text entry of the isolation
           sources as well as on the pubmed identifiers from which the GI obtained.

    * `relative`: Should we divide the counts of every input sequence by the
                      number of envo terms that were associated to it.
                      Defaults to `True`.

    * `proportional`: Should we divide the counts of every input sequence by the
                      number of envo terms that were associated to it.
                      Defaults to `True`.

    * `backtracking`: For every term identified by the tagger, we will propagate
                      frequency counts up the acyclic directed graph described by
                      the ontology. Defaults to `False`.

    * `restrict`: Restrict the output to the descendants of just one ENVO term.
                  This removes all other terms that are not reachable through the node.
                  For instance you could specify: `ENVO:00010483`. Disabled by default.

    * `num_threads`: The number of threads. Default to the number of cores on the
                     current machine or a maximum of 32.

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
                 normalization = 'flat',
                 proportional  = True,
                 backtracking  = False,
                 restrict      = None,
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
        # Restrict parameter #
        message = "The '--restrict' parameter must be an ENVO term, not '%s'."
        if restrict and not restrict[:5] == 'ENVO:': raise Exception(message % restrict)
        message = "The '--restrict' parameter must be a known ENVO term."
        if restrict and not restrict in self.serial_to_concept.values(): raise Exception(message)
        self.restrict = restrict
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
        if git_repo: print "The exact version of the codebase is: " + git_repo.short_hash
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

    @property_cached
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

    @property_cached
    def only_top_sequences(self):
        """Make a new fasta file where only the top N sequences are included
        (in terms of their abundance). Skipped if no abundance info is given."""
        if not self.abundances: return self.renamed_fasta
        if self.N is None:      return self.renamed_fasta
        # Parse it #
        N = int(self.N)
        # Create file #
        only_top_fasta = FASTA(self.out_dir + 'top_seqs.fasta')
        # Print status #
        print "Using: " + self.renamed_fasta
        print "--> STEP 1B: Get the top %i sequences (in terms of their abundances)." % N
        # Check the user inputted value #
        if N > self.input_file.count:
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
    # gis_with_envos
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
    def db(self):
        """The sqlite3 database containing every GI number in all NCBI that has an
        isolation source associated to it. In addition, the pubmed-ID is listed too
        if there is one. If we don't have it locally already, we will go download it.
        Thus, the database containing two tables:

        CREATE TABLE "gi"
            (
              "id"      INTEGER PRIMARY KEY NOT NULL,
              "isokey"  INTEGER NOT NULL REFERENCES "isolation"("id"),
              "pubmed"  INTEGER
            );
        CREATE TABLE "isolation"
            (
              "id"     INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
              "source" TEXT NOT NULL,
              "envos"  BLOB NOT NULL <- marshaled tuple of ints
            );
        CREATE VIEW "data" AS SELECT gi.id, (->source), gi.pubmed, (->envos) from "gi";
        CREATE INDEX "gi_index" on "gi" (id);
        CREATE INDEX "isolation_index" on "isolation" (source);
        """
        print "--> STEP 5: Setting up the SQLite3 database connection."
        path     = module_dir + 'data_sources/gi_db.sqlite3'
        drop_box = "1yxoc1530am6tn0/gi_db.sqlite3.zip"
        retrieve = "https://www.dropbox.com/s/%s?dl=1" % drop_box
        md5      = "bda6297a8cc6219836ed1cb118e03510"
        database = Database(path, retrieve=retrieve, md5=md5, text_fact=bytes)
        self.timer.print_elapsed()
        return database

    @property_cached
    def gis_with_envo(self):
        """A set containing every GI that was found and that had some envo numbers associated."""
        return set(gi for gi in self.unique_gis if gi in self.db)

    @property_cached
    def seq_to_counts(self):
        """A dictionary linking every input sequence to its summed normalized concept
        counts dict, provided the input sequence had some hits, and at least one hit had
        some envo number associated. Otherwise it is absent.
        NB: What we want to account for is the fact that two GIs originating from the
        same sequence could be pointing to the same isolation source.
        In such case, we shall count the envo numbers from that isolation source only once.
        We can also avoid counting two GIs that are coming from the same study, if a pubmed
        number is available."""
        # Message #
        msg = "Got %i GI hits and %i of them had one or more EnvO terms associated."
        print msg % (len(self.unique_gis), len(self.gis_with_envo))
        # Step 5 #
        print "--> STEP 6: Computing EnvO term frequencies."
        results = {}
        # Useful later #
        gi_to_key    = lambda gi: self.db.get("gi","id",gi)[1]
        gi_to_pub    = lambda gi: self.db.get("gi","id",gi)[2]
        key_to_envos = lambda k: marshal.loads(self.db.get("isolation","id",k)[2])
        first_key    = lambda d,k: [(gi, knp) for gi, knp in d.items()  if knp[0]==k][0]
        first_pub    = lambda d,p: [(gi, knp) for gi, knp in d.items() if knp[1]==p][0]
        # "Flat" or "Unique source" are quite similar #
        set_or_list = None
        if self.normalization == 'flat': set_or_list = list
        if self.normalization == 'ui':   set_or_list = set
        # It's just about using a list or a set in the right place #
        if self.normalization == 'flat' or self.normalization == 'ui':
            for seq, gis in self.seq_to_gis.items():
                gis_w_envo = [gi for gi in gis if gi in self.db]
                if not gis_w_envo: continue
                isokeys   = set_or_list(gi_to_key(gi) for gi in gis_w_envo)
                all_envos = [e for key in isokeys for e in key_to_envos(key)]
                if self.backtracking:
                    all_envos.extend([p for e in all_envos for p in self.child_to_parents[e]])
                if self.restrict:
                    all_envos = [e for e in all_envos if self.ontology.descends(e, self.restrict)]
                if not all_envos: continue
                if self.proportional: score = 1/len(all_envos)
                else:                 score = 1.0
                counts = defaultdict(float)
                for e in all_envos: counts[e] += score
                results[seq] = counts
        # Unique source and unique pubmed #
        if self.normalization == 'upui':
            for seq, gis in self.seq_to_gis.items():
                gis_w_envo = [gi for gi in gis if gi in self.db]
                if not gis_w_envo: continue
                gi_to_knp  = {gi: (gi_to_key(gi), gi_to_pub(gi)) for gi in gis_w_envo}
                uniq_keys  = set(key for key, pub in gi_to_knp.values())
                gi2knp_uk  = dict(first_key(gi_to_knp, key) for key in uniq_keys)
                pubmeds    = set(pub for key, pub in gi2knp_uk.values() if pub != None)
                gi2knp_up  = dict(first_pub(gi2knp_uk, pubmed) for pubmed in pubmeds)
                # But if the pubmed field is None, keep them all #
                gi2knp_up.update({(gi, knp) for gi, knp in gi2knp_uk.items() if knp[1]==None})
                all_envos  = [e for key,pub in gi2knp_up.values() for e in key_to_envos(key)]
                if self.backtracking:
                    all_envos.extend([p for e in all_envos for p in self.child_to_parents[e]])
                if self.restrict:
                    all_envos = [e for e in all_envos if self.ontology.descends(e, self.restrict)]
                if not all_envos: continue
                if self.proportional: score = 1/len(all_envos)
                else:                 score = 1.0
                counts = defaultdict(float)
                for e in all_envos: counts[e] += score
                results[seq] = counts
        # Check #
        if not results: raise Exception("We found no isolation sources with your input. Sorry.")
        # Return #
        self.timer.print_elapsed()
        return results

    # --------------------------- In this section --------------------------- #
    # serial_to_concept
    # child_to_parents
    # concept_to_name
    # ontology

    @class_property
    @classmethod
    @cached
    def serial_to_concept(cls):
        """A dictionary linking every concept serial to its concept id.
        Every line in the file contains three columns: serial, concept_type, concept_id.
        An example line: 1007000022 -27 ENVO:00000021"""
        return dict(line.split()[0::2] for line in open(module_dir + 'data_envo/envo_entities.tsv'))

    @class_property
    @classmethod
    @cached
    def child_to_parents(cls):
        """A dictionary linking concept id to a list of parent concept ids.
        Every line in the file contains two columns: child_serial, parent_serial
        An example line: 1007000003 1007001640"""
        result = defaultdict(list)
        with open(module_dir + 'data_envo/envo_groups.tsv') as handle:
            for line in handle:
                child_serial, parent_serial = line.split()
                child_concept = cls.serial_to_concept[child_serial]
                parent_concept = cls.serial_to_concept[parent_serial]
                result[child_concept].append(parent_concept)
        return result

    @class_property
    @classmethod
    @cached
    def concept_to_name(cls):
        """A dictionary linking the concept id to relevant names. In this case ENVO terms.
        Hence, ENVO:00000095 would be linked to 'lava field'
        An example line: ENVO:00000015  ocean"""
        handle = open(module_dir + 'data_envo/envo_preferred.tsv')
        result = {}
        for line in handle:
            envo, name = line.strip('\n').split('\t')
            if envo == "ENVO:root": continue
            result[envo] = name
        return result

    @class_property
    @classmethod
    @cached
    def integer_to_name(cls):
        """A dictionary linking the concept id integer to relevant names.
        The same thing as above, but with integers as keys."""
        handle = open(module_dir + 'data_envo/envo_preferred.tsv')
        result = {}
        for line in handle:
            envo, name = line.strip('\n').split('\t')
            if envo == "ENVO:root": continue
            result[int(envo[5:])] = name
        return result

    @class_property
    @classmethod
    @cached
    def ontology(self):
        """The ontology instance (singleton).
        Give access to all the ENVO terms graph."""
        return Ontology()
