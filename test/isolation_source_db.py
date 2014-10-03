"""
=================================
Test the NCBI isolation source db
=================================
"""

################################################################################
from seqenv import Analysis
d = "examples/samples/"
a = Analysis(d+'community.fasta', out_dir=d+"output/", abundances=d+"abundances.tsv", N=10)

via_db = {gi: a.gi_to_text.get(gi) for gi in a.gi_to_matches.keys()}

from eutils import gis_to_records, record_to_source
recs = r = gis_to_records(a.gi_to_matches.keys())
via_eutils = {k:record_to_source(v) for k,v in recs.items()}

print set(via_db) == set(via_eutils)
all([via_db.get(k) == via_eutils.get(k) for k in via_db])