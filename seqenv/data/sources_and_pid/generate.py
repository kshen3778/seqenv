#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
A script to parse and download information from the NCBI NT database
useful for the seqenv project.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ ./generate.py
"""

# Built-in modules #
import os, time, sqlite3
from shell_command import shell_output

# Internal modules #
from seqenv.common.timer import Timer

# Third party modules #
from Bio import Entrez
from Bio.Entrez.Parser import CorruptedXMLError
from tqdm import tqdm

# Constants #
Entrez.email = "I don't know who will be running this script"
at_a_time = 40

################################################################################
def gis_to_records(gis, progress=False):
    """Download information from NCBI in batch mode
    Input is a list of Gids."""
    # Should we display progress ? #
    progress = tqdm if progress else lambda x:x
    # Do it by chunks #
    gis = list(gis)
    result = {}
    # Main loop #
    for i in progress(range(0, len(gis), at_a_time)):
        chunk = gis[i:i+at_a_time]
        records = chunk_to_records(chunk)
        # Only store the things we need in memory, otherwise we will get into swap space #
        sources = map(record_to_source, records)
        pubmed_ids =  map(record_to_pubmed_id, records)
        # Update dictionary #
        result.update(dict(zip(chunk, zip(sources, pubmed_ids))))
    # Return #
    return result

def chunk_to_records(chunk):
    """Download from NCBI until it works. Will restart until reaching the python
    recursion limit. We don't want to get our IP banned from NCBI so we have
    a little pause at every function call."""
    time.sleep(0.5)
    try:
        response = Entrez.efetch(db="nucleotide", id=chunk, retmode="xml")
        records = list(Entrez.parse(response, validate=True))
        return records
    except CorruptedXMLError:
        return chunk_to_records(chunk)

def record_to_source(record):
    qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in qualifiers:
        if qualifier['GBQualifier_name'] == 'isolation_source':
            return qualifier['GBQualifier_value']

def record_to_pubmed_id(record):
    references = record['GBSeq_references'][0]
    pubmed_id = references.get('GBReference_pubmed')
    return pubmed_id

###############################################################################
def result_dict_to_lines(results):
    for gi, info in results.iteritems():
        yield gi + ',' + str(info[0]) + ',' + str(info[1]) + '\n'

def serialize_results(results, sqlite_file):
    if os.path.exists(sqlite_file): os.remove(sqlite_file)
    with sqlite3.connect(sqlite_file) as connection:
        cursor = connection.cursor()
        cursor.execute("CREATE table 'data' (gi integer, source text, pubid integer)")
        sql_command = "INSERT into 'data' values (?,?,?)"
        values = ((gi, info[0], info[1]) for gi, info in results.iteritems())
        cursor.executemany(sql_command, values)
        cursor.execute("CREATE INDEX if not exists 'data_index' on 'data' (gi)")
        connection.commit()
        cursor.close()

def test():
    """Test this script"""
    # Just a bunch GIs numbers #
    test_gis = ['6451693', '127', '76365841', '22506766', '389043336',
                '497', '429143984', '264670502', '74268401', '324498487']
    # The result we should get #
    template_result = """
    74268401,None,None
    6451693,None,None
    76365841,Everglades wetlands,16907754
    324498487,bacterioplankton sample from lake,None
    389043336,lake water at 5 m depth during dry season,None
    22506766,None,None
    127,None,None
    429143984,downstream along river bank,None
    497,None,3120795
    264670502,aphotic layer; anoxic zone; tucurui hydroeletric power plant reservoir,None
    """
    # Make it pretty #
    template_result = '\n'.join(l.lstrip(' ') for l in template_result.split('\n') if l)
    # The result we got #
    results = gis_to_records(test_gis)
    # Check #
    text_version = ''.join(result_dict_to_lines(results))
    assert text_version == template_result
    # Return #
    return results

###############################################################################
if __name__ == '__main__':
    # Get current directory #
    current_dir = os.path.dirname(os.path.realpath(__file__)) + '/'

    # Start timer #
    timer = Timer()
    timer.print_start()

    # Test the pipeline #
    print 'STEP 0: Testing the script and connection'
    test()

    # Get all the G-ids from the current local NT database #
    print 'STEP 1: Get all GIs from local nt database into a file (about 2h)'
    gids_file = current_dir + 'all_gis.txt'
    shell_output("blastdbcmd -db nt -entry all -outfmt '%g' > " + gids_file)
    timer.print_elapsed()

    # Parse the file we created #
    print 'STEP 2: Parse file created'
    with open(gids_file, 'r') as handle: all_gis = list(handle)
    timer.print_elapsed()

    # Now we query NCBI with the eutils to get every isolation source and pubmed ID #
    # This will run for very long #
    print 'STEP 3: Query NCBI for all records (about 1 or 2 days)'
    results = gis_to_records(all_gis, progress=tqdm)
    timer.print_elapsed()

    # Serialize the result to disk in an sqlite databse #
    print 'STEP 4: Writing results from memory to disk in an sqlite databse'
    sqlite_file = current_dir + 'gi_to_source.tsv'
    serialize_results(results, sqlite_file)
    timer.print_elapsed()

    # Clean up #
    #os.remove(gids_file)

    # End #
    print 'Done results are in "%s"' % os.path.abspath(current_dir)
    timer.print_end()
    timer.print_total_elapsed()