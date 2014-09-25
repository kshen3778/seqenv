"""
========================
Generate fake abundances
========================
"""

# Modules #
import os, inspect, numpy, pandas, names
from seqenv.fasta import FASTA

# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = FASTA(current_dir + "community.fasta")

################################################################################
def data():
    """Create some fake data in a dataframe"""
    x_size = len(fasta)
    y_size = 10
    numpy.random.seed(0)
    M = numpy.random.randint(0, 1000, (x_size, y_size))
    df = pandas.DataFrame(M, index=[seq.id for seq in fasta], columns=[names.get_first_name() for j in range(y_size)])
    return df


df = data()
df = df.apply(lambda x: x/x.sum())
df.to_csv(current_dir + 'abundances.tsv', sep='\t')