"""
=========================
Run the sequences example
=========================
"""

# Modules #
from seqenv import Analysis

################################################################################
analysis = Analysis("community.fasta", out_dir='output/', num_threads=20)

# Run #
analysis.timer.print_start()
analysis.outputs.make_all()
analysis.timer.print_end()
analysis.timer.print_total_elapsed()