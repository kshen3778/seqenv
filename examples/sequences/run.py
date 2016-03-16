"""
=========================
Run the sequences example
=========================
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
analysis = Analysis("community.fasta", out_dir='output/', num_threads=20)

# Run #
analysis.run()