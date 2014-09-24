"""
=======================
Test the whole pipeline
=======================
"""

# Modules #
import os, inspect
from seqenv import Analysis

# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = current_dir + "community.fasta"

################################################################################
analysis = Analysis(fasta, num_threads=3, min_identity=0.1, e_value=0.1, backtracking=True)
#analysis.gi_to_concepts