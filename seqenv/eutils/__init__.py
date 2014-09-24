# Built-in modules #

# Third party modules #
from Bio import Entrez
Entrez.email = "I don't know who will be running this script"

################################################################################
def gi_to_source(gi):
    gb_entry = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="xml")
    gb_record = Entrez.read(gb_entry)[0]
    gb_qualifiers = gb_record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in gb_qualifiers:
        key, value = qualifier.items()
        if key == 'isolation_source': return value

################################################################################
def gi_to_abstract(gi):
    gb_entry = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="xml")
    gb_record = Entrez.read(gb_entry)[0]
    gb_qualifiers = gb_record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in gb_qualifiers:
        key, value = qualifier.items()
        if key == 'isolation_source': return value