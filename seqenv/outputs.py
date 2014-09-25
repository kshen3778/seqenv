# Built-in modules #

# Internal modules #
from seqenv.common.cache import property_cached

# Third party modules #
import pandas

################################################################################
class OutputGenerator(object):
    """Once the Analysis is run and all the data is in memory in the
    form of python objects, this class will take care of generating
    all the output files the user could possibly want. You pass it the Analysis
    object obviously."""

    sep = '\t'
    float_format = '%.5g'

    def __init__(self, analysis):
        self.analysis = analysis
        self.out_dir = analysis.out_dir

    def make_all(self):
        """Let's generate all the files"""
        with open(self.out_dir + 'seq_to_concepts.csv', 'w') as handle: handle.writelines(self.csv_seq_to_concepts())
        with open(self.out_dir + 'seq_to_names.csv', 'w')    as handle: handle.writelines(self.csv_seq_to_names())

    @property_cached
    def df_seqs_concepts(self):
        """A normalized matrix with sequences as columns and concepts as rows."""
        # Get the data #
        df = pandas.DataFrame(self.analysis.seq_to_counts)
        df = df.fillna(0)
        # Rename to original names #
        renamed_to_orig = dict((v,k) for k,v in self.analysis.orig_names_to_renamed.iteritems())
        df = df.rename(columns=renamed_to_orig)
        # Return
        return df

    def csv_seq_to_concepts(self):
        """A CSV file"""
        return self.df_seqs_concepts.to_csv(None, sep=self.sep, float_format=self.float_format)

    def csv_seq_to_names(self, sep='\t'):
        """A CSV file"""
        df = self.df_seqs_concepts.rename(index=self.analysis.concept_to_name)
        return df.to_csv(sep=self.sep, float_format=self.float_format)

    def csv_otu_to_names(self, sep='\t'):
        """A CSV file"""
        # Get both matrices #
        df1 = self.df_seqs_concepts.rename(index=self.analysis.concept_to_name).transpose()
        df2 = self.abundance_df.transpose()
        # Multiply them #
        df = df1.dot(df2)
        return df.to_csv(sep=self.sep, float_format=self.float_format)

    def output_1(self):
        """The counts per sequence, one concept per line:
        - OTU1, ENVO:00001, 4, GIs : [56, 123, 345]
        - OTU1, ENVO:00002, 7, GIs : [56, 123, 345]
          or
        - OTU1, ocean, 4, GIs : [56, 123, 345]
        """
        for seq, counts in self.seq_to_counts.items():
            '\t'.join(seq)

    def output_2(self):
        """Possible wanted output 2:
        - GI56 : oceanic soil
          or
        - GI56 : PubMed6788
        """
        pass

    def output_3(self):
        """ Possible output 3:
        - ENVO:00001 : ocean
        """
        pass
