"""
==================
Test the SeqSearch
==================
"""

# Modules #
import os, inspect
from seqenv.fasta import FASTA
from seqenv.seqsearch.parallel import ParallelSeqSearch

# Constants #
current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
fasta = FASTA(current_dir + "../examples/samples/community.fasta")

################################################################################
search = ParallelSeqSearch(fasta, 'nucl', 'nt', algorithm='blast', num_threads=3)
search.run()