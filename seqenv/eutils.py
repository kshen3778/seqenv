# Built-in modules #

# Third party modules #
from Bio import Entrez
from tqdm import tqdm

# Constants #
Entrez.email = "I don't know who will be running this script"

################################################################################
def gis_to_records(gis, progress=True):
    # Should we display progress ? #
    progress = tqdm if progress else lambda x:x
    # Do it by chunks #
    gis = list(gis)
    at_a_time = 1000
    result = {}
    # Main loop #
    for i in progress(range(0, len(gis), at_a_time)):
        chunk = gis[i:i+at_a_time]
        entries = Entrez.efetch(db="nucleotide", id=chunk, rettype="gb", retmode="xml")
        result.update(dict(zip(gis, Entrez.parse(entries, validate=True))))
    # Return #
    return result

################################################################################
def record_to_source(record):
    qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in qualifiers:
        if qualifier['GBQualifier_name'] == 'isolation_source':
            return qualifier['GBQualifier_value']

################################################################################
def record_to_abstract(record):
    raise NotImplemented('')