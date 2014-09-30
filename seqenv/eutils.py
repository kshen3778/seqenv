# Built-in modules #

# Third party modules #
import time
from Bio import Entrez
from Bio.Entrez.Parser import CorruptedXMLError
from tqdm import tqdm

# Constants #
Entrez.email = "I don't know who will be running this script"

################################################################################
def gis_to_records(gis, progress=False):
    """Download information from NCBI in batch mode"""
    # Should we display progress ? #
    progress = tqdm if progress else lambda x:x
    # Do it by chunks #
    gis = list(gis)
    at_a_time = 40
    result = {}
    # Main loop #
    for i in progress(range(0, len(gis), at_a_time)):
        chunk = gis[i:i+at_a_time]
        records = chunk_to_records(chunk)
        result.update(dict(zip(chunk, records)))
    # Return #
    return result

################################################################################
def chunk_to_records(chunk):
    """Download from NCBI until it works. Will restart until reaching the python
    recursion limit. We don't want to get banned from NCBI so we have a little
    pause at every function call."""
    time.sleep(0.5)
    print chunk
    try:
        print "-"
        response = Entrez.efetch(db="nucleotide", id=chunk, retmode="xml")
        print "="
        records = list(Entrez.parse(response, validate=True))
        return records
    except CorruptedXMLError:
        print "*"
        return chunk_to_records(chunk)

################################################################################
def record_to_source(record):
    qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in qualifiers:
        if qualifier['GBQualifier_name'] == 'isolation_source':
            return qualifier['GBQualifier_value']

################################################################################
def record_to_abstract(record):
    raise NotImplemented('')