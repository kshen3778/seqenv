# Built-in modules #
import warnings, marshal, re

# Internal modules #
import seqenv
from seqenv.common.cache     import property_cached
from seqenv.common.autopaths import DirectoryPath

# Third party modules #
import pandas

# We don't want the annoying h5py warning #
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import biom

################################################################################
class OutputGenerator(object):
    """Once the Analysis is done running and all the data is in memory in the
    form of python objects, this object will take care of generating
    all the output files the user could possibly want. You pass it the Analysis
    object obviously."""

    sep = '\t'
    float_format = '%.5g'

    def __init__(self, analysis):
        self.analysis = analysis
        self.a        = analysis

    # --------------------------- In this section --------------------------- #
    # make_all
    # df_seqs_concepts
    # abundance_mat_mult
    # df_sample_concepts
    # df_sample_names

    def make_all(self):
        """Let's generate all the files"""
        # General matrices #
        self.tsv_seq_to_concepts()
        self.tsv_seq_to_names()
        self.list_sequence_concept()
        # Only in the with 'samples' case #
        if self.a.abundances: self.tsv_samples_to_names()
        if self.a.abundances: self.biom_output()
        # Graphical outputs #
        self.per_seq_dot_files()
        if self.a.abundances: self.per_sample_dot_files()

    @property_cached
    def df_seqs_concepts(self):
        """A normalized dataframe with sequences as columns and concepts (envo terms) as rows."""
        # Get the data #
        df = pandas.DataFrame(self.a.seq_to_counts)
        df = df.fillna(0)
        # Rename to original names #
        df = df.rename(columns=self.a.renamed_to_orig)
        # Rename envo integers to envo strings #
        envo_int_to_id = lambda e: "ENVO:%08d" % e
        df = df.rename(index=envo_int_to_id)
        # Return
        return df

    def abundance_mat_mult(self, human_names=False):
        """We operate a matrix multiplication with the abundances
        file provided, to link samples to either concept human readable names
        or the envo term IDs."""
        # Get results #
        df1 = self.df_seqs_concepts
        # Rename #
        if human_names: df1 = df1.rename(index=self.a.concept_to_name)
        # Remove those that were discarded #
        df2 = self.a.df_abundances
        df2 = df2.loc[df1.columns]
        # Odd bug detection #
        if any(s[:1].isdigit() for s in df2.columns):
            msg = "None of the sample names in file '%s' can start with a number"
            raise Exception(msg % self.a.abundances.filename)
        # Multiply them (dot product) #
        assert all(df1.columns == df2.index)
        df = df1.dot(df2)
        # Return
        return df

    @property_cached
    def df_sample_concepts(self):
        """A dataframe where we operate a matrix multiplication with the abundances
        file provided, to link samples to concept envo terms."""
        return self.abundance_mat_mult()

    @property_cached
    def df_sample_names(self):
        """A dataframe where we operate a matrix multiplication with the abundances
        file provided, to link samples to concept human readable names."""
        return self.abundance_mat_mult(True)

    # --------------------------- In this section --------------------------- #
    # tsv_seq_to_concepts
    # tsv_seq_to_names
    # tsv_samples_to_names
    # list_sequence_concept
    # biom_output

    def tsv_seq_to_concepts(self, name="seq_to_concepts.tsv"):
        """A TSV matrix file containing the df_seqs_concepts matrix"""
        with open(self.a.out_dir + name, 'w') as handle:
            content = self.df_seqs_concepts.to_csv(None, sep=self.sep, float_format=self.float_format)
            handle.writelines(content)

    def tsv_seq_to_names(self, name='seq_to_names.tsv'):
        """A TSV matrix file where we translate the concept to human readable names"""
        with open(self.a.out_dir + name, 'w') as handle:
            df = self.df_seqs_concepts.rename(index=self.a.concept_to_name)
            content = df.to_csv(None, sep=self.sep, float_format=self.float_format)
            handle.writelines(content)

    def tsv_samples_to_names(self, name='samples_to_names.tsv'):
        """A TSV matrix file with matrix `df_sample_names`."""
        with open(self.a.out_dir + name, 'w') as handle:
            content = self.df_sample_names.to_csv(None, sep=self.sep, float_format=self.float_format)
            handle.writelines(content)

    def list_sequence_concept(self, name='list_concepts_found.tsv'):
        """A flat TSV file listing every concept found for every sequence.
        It has one concept per line and looks something like this:
        - OTU1, ENVO:00001, ocean, 4, GIs : [56, 123, 345]
        - OTU1, ENVO:00002, soil, 7, GIs : [22, 44]
        """
        # Useful later #
        gi_to_key    = lambda gi: self.a.db.get("gi","id",gi)[1]
        key_to_envos = lambda key: marshal.loads(self.a.db.get("isolation","id",key)[2])
        gi_to_envos  = lambda gi: key_to_envos(gi_to_key(gi))
        # Loop #
        with open(self.a.out_dir + name, 'w') as handle:
            for seq, gis in self.a.seq_to_gis.items():
                gis     = [gi for gi in gis if gi in self.a.db]
                isokeys = set(gi_to_key(gi) for gi in gis)
                envos   = [e for key in isokeys for e in key_to_envos(key)]
                for envo in envos:
                    seq_name     = self.a.renamed_to_orig[seq]
                    envo_id      = "ENVO:%08d" % envo
                    concept_name = self.a.integer_to_name.get(envo, envo_id)
                    concept_gis  = [gi for gi in gis if envo in gi_to_envos(gi)]
                    count_gis    = len(concept_gis)
                    line         = (seq_name, envo_id, concept_name, str(count_gis), str(concept_gis))
                    handle.write('\t'.join(line) + '\n')

    def biom_output(self, name='samples.biom'):
        """The same matrix as the user gave in the abundance file, but with source
        information attached for every sequence.
        See http://biom-format.org"""
        data = self.a.df_abundances
        with open(self.a.out_dir + name, 'w') as handle:
            # Basic #
            sample_ids = data.columns
            sample_md = None
            observation_ids = data.index
            # Observation metadata #
            observation_md = []
            for seq in data.index:
                seq_name = self.a.orig_names_to_renamed[seq]
                counts = self.a.seq_to_counts.get(seq_name)
                if not counts: observation_md.append({})
                else: observation_md.append({'source': counts})
            # Output #
            t = biom.table.Table(data.transpose().as_matrix(), sample_ids, observation_ids, sample_md, observation_md)
            handle.write(t.to_json('seqenv version %s') % seqenv.__version__)

    # --------------------------- In this section --------------------------- #
    # per_seq_dot_files
    # per_sample_dot_files

    def per_seq_dot_files(self):
        """Generations of files that can be viewed in `graphviz`.
        There is one dotfile per every input sequence.
        We also automiatcally make a corresponding PDF file."""
        # The output directory #
        directory = DirectoryPath(self.a.out_dir+'per_seq_ontology/')
        directory.create_if_not_exists()
        # Main loop #
        for seq in self.a.seq_to_counts:
            dot_path = directory + seq + '.dot'
            pdf_path = directory + seq + '.pdf'
            counts = self.a.seq_to_counts[seq]
            counts = {"ENVO:%08d"%k:v for k,v in counts.items()}
            total  = sum(counts.values())
            counts = {k:v/total for k,v in counts.items()}
            envos  = counts.keys()
            graph  = self.a.ontology.get_subgraph(envos)
            graph  = self.a.ontology.add_weights(graph, counts)
            graph  = self.a.ontology.add_style(graph)
            self.a.ontology.write_to_dot(graph, dot_path)
            self.a.ontology.add_legend(dot_path)
            self.a.ontology.draw_to_pdf(dot_path, pdf_path)

    def per_sample_dot_files(self):
        """Generations of files that can be viewed in `graphviz`.
        There is one dotfile per every sample inputted.
        We also automiatcally make a corresponding PDF file."""
        # The output directory #
        directory = DirectoryPath(self.a.out_dir+'per_sample_ontology/')
        directory.create_if_not_exists()
        # Main loop #
        for i, sample in self.df_sample_concepts.iteritems():
            # File path #
            sanitized_name = "".join([c for c in sample.name if re.match(r'\w', c)])
            dot_path = directory + sanitized_name +'.dot'
            pdf_path = directory + sanitized_name +'.pdf'
            # Counts #
            counts = sample / sample.sum()
            counts = dict(counts)
            envos  = counts.keys()
            # Skip sample if it has no counts #
            if sample.sum() == 0: continue
            # Make graph #
            graph  = self.a.ontology.get_subgraph(envos)
            graph  = self.a.ontology.add_weights(graph, counts)
            graph  = self.a.ontology.add_style(graph)
            # Write output #
            self.a.ontology.write_to_dot(graph, dot_path)
            self.a.ontology.add_legend(dot_path)
            self.a.ontology.draw_to_pdf(dot_path, pdf_path)

    # --------------------------- In this section --------------------------- #
    # output_3

    def output_3(self):
        """Possible output #3: the number of terms per OTU
        OTU1: 0
        OTU2: 2
        OTU3: 1"""
        pass