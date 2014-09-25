"""
========================
Generate fake abundances
========================
"""

# Modules #
import os, inspect, numpy, random, scipy, pandas, names
from seqenv.fasta import FASTA

# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = FASTA(current_dir + "community.fasta")

################################################################################
def data():
    """Create some fake data in a dataframe"""
    numpy.random.seed(0)
    random.seed(0)
    x_size = len(fasta)
    y_size = 10
    data = scipy.rand(max(x_size,y_size))
    M = scipy.zeros([x_size,y_size])
    for i in range(x_size):
        for j in range(y_size): M[i,j] = abs(data[i] - data[j])
    df = pandas.DataFrame(M, index=[seq.id for seq in fasta], columns=[names.get_first_name() for j in range(y_size)])
    return df

df = data()
df.to_csv(current_dir + 'abundances.tsv', sep='\t')