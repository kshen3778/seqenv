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
out_dir = current_dir + 'output/'

################################################################################
analysis = Analysis(fasta, out_dir=out_dir)
analysis.timer.print_start()
analysis.outputs.make_all()
analysis.timer.print_end()
analysis.timer.print_total_elapsed()