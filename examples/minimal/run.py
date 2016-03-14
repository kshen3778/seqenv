"""
=======================
Run the minimal example
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
analysis = Analysis("test.fasta", out_dir='output/', num_threads=1)

# Run #
analysis.run()