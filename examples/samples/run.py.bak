"""
=======================
Run the samples example
=======================
"""

# Modules #
from seqenv import Analysis

# Change directory #
import inspect, os
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'
os.chdir(current_dir)

################################################################################
# Main object #
analysis = Analysis("community.fasta", out_dir='output/',
                    abundances="abundances.tsv", N=10, num_threads=5)

# Run #
analysis.run()