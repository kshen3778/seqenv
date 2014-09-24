"""
===============================
Test the command line interface
===============================
"""

# Modules #
import os, inspect, sh

# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = current_dir + "community.fasta"

################################################################################
sh.seqenv(fasta)