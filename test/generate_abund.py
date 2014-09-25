"""
========================
Generate fake abundances
========================
"""

# Modules #
import os, inspect, numpy, random, scipy, pandas, names


# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = current_dir + "community.fasta"

################################################################################
def data(self):
    """Create some fake data in a dataframe"""
    numpy.random.seed(0)
    random.seed(0)
    x_size =
    y_size =
    x = scipy.rand(size)
    M = scipy.zeros([size,size])
    for i in range(size):
        for j in range(size): M[i,j] = abs(x[i] - x[j])
    df = pandas.DataFrame(M, index=[names.get_last_name() for _ in range(size)],
                             columns=[names.get_first_name() for _ in range(size)])
    df['Mary']['Day'] = 1.5
    df['Issac']['Day'] = 1.0
    return df

df = data()
df.to_csv(current_dir + 'abundances.tsv', sep='\t')