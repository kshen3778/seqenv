#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
A script to parse and download information from the NCBI NT database
useful for the seqenv project. We want isolation sources and pubmed_ids
for all GI numbers found in the NT database.

See https://biosupport.se/p/719/

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ ./generate.py
"""

# Built-in modules #
import os, time, sqlite3, inspect
from itertools import islice
import urllib2

# Internal modules #
from seqenv.common.timer import Timer
from seqenv.common.autopaths import FilePath

# Third party modules #
from Bio import Entrez
from Bio.Entrez.Parser import CorruptedXMLError
from tqdm import tqdm
from shell_command import shell_output

# Constants #
Entrez.email = "I don't know who will be running this script"
at_a_time = 100

# Get current directory #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename))  + '/'

# Files #
gids_file   = FilePath(current_dir + 'all_gis.txt')
sqlite_file = FilePath(current_dir + 'gi_index.sqlite3')

###############################################################################
def all_gis():
    # If the gi files doesn't exist #
    if gids_file.count_bytes < 1:
        print 'STEP 1: Get all GIs from local nt database into a file (about 6h).'
        shell_output("blastdbcmd -db nt -entry all -outfmt '%g' > " + gids_file)
    else:
        print '-> All GI numbers already found in file "%s", skipping STEP 1.' % gids_file.filename
    print "-> Got %i total GI numbers" % len(gids_file)
    # If the database doesn't exists #
    if sqlite_file.count_bytes < 1:
        print "-> The database doesn't exist. Creating it."
        return gids_file.lines, len(gids_file)
    else:
        print '-> The database seems to exist already. Attempting to restart where it was stopped.'
        connection = sqlite3.connect(sqlite_file)
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM data ORDER BY ROWID DESC LIMIT 1;")
        last_gi = str(cursor.next()[0])
        cursor.close()
        connection.close()
        generator = gids_file.lines
        n = 0
        while True:
            gi = generator.next()
            n += 1
            if gi == last_gi: break
        else:
            raise Exception("Didn't find the last GI in the all_gis.txt file")
        return generator, len(gids_file) - n

###############################################################################
def gis_to_records(gids, verbose=True):
    """Download information from NCBI in batch mode"""
    if verbose: print 'STEP 2: Querying NCBI and writing to database (about 300h)'
    progress = tqdm if verbose else lambda x:x
    for i in progress(xrange(0, len(gids_file), at_a_time)):
        chunk      = list(islice(gids, 0, at_a_time))
        records    = chunk_to_records(chunk)
        sources    = map(record_to_source, records)
        pubmed_ids = map(record_to_pubmed_id, records)
        # Generate the result #
        has_source = [i for i in xrange(len(chunk)) if sources[i]]
        if not has_source: continue
        yield {chunk[i]: (sources[i], pubmed_ids[i]) for i in has_source}

def chunk_to_records(chunk):
    """Download from NCBI until it works. Will restart until reaching the python
    recursion limit. We don't want to get our IP banned from NCBI so we have
    a little pause at every function call."""
    time.sleep(0.2)
    try:
        response = Entrez.efetch(db="nucleotide", id=chunk, retmode="xml")
        records = list(Entrez.parse(response, validate=True))
        return records
    except (CorruptedXMLError, urllib2.HTTPError, urllib2.URLError):
        time.sleep(5)
        return chunk_to_records(chunk)

def record_to_source(record):
    qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
    for qualifier in qualifiers:
        if qualifier['GBQualifier_name'] == 'isolation_source':
            return qualifier['GBQualifier_value']

def record_to_pubmed_id(record):
    if 'GBSeq_references' not in record: return None
    references = record['GBSeq_references'][0]
    pubmed_id = references.get('GBReference_pubmed')
    return pubmed_id

###############################################################################
def add_to_database(results):
    connection = sqlite3.connect(sqlite_file, isolation_level=None)
    cursor = connection.cursor()
    cursor.execute("CREATE table 'data' (gi integer, source text, pubid integer)")
    sql_command = "INSERT into 'data' values (?,?,?)"
    for chunk in results:
        values = ((gi, info[0], info[1]) for gi, info in chunk.iteritems())
        cursor.executemany(sql_command, values)
    cursor.execute("CREATE INDEX if not exists 'data_index' on 'data' (gi)")
    connection.commit()
    cursor.close()
    connection.close()

###############################################################################
def test():
    """Test this script"""
    # Just a bunch of random GI numbers #
    test_gis = ['6451693', '127', '76365841', '22506766', '389043336',
                '497', '429143984', '264670502', '74268401', '324498487']
    # The result we should get #
    template_result = """
    429143984,downstream along river bank,None
    76365841,Everglades wetlands,16907754
    324498487,bacterioplankton sample from lake,None
    389043336,lake water at 5 m depth during dry season,None
    264670502,aphotic layer; anoxic zone; tucurui hydroeletric power plant reservoir,None
    """
    # Make it pretty #
    template_result = '\n'.join(l.lstrip(' ') for l in template_result.split('\n') if l)
    # The result we got #
    results = gis_to_records(iter(test_gis), verbose=False).next()
    # Function #
    def result_dict_to_lines(results):
        for gi, info in results.items():
            yield gi + ',' + str(info[0]) + ',' + str(info[1]) + '\n'
    # Check #
    text_version = ''.join(result_dict_to_lines(results))
    assert text_version == template_result
    # Return #
    return results

###############################################################################
if __name__ == '__main__':
    # Start timer #
    timer = Timer()
    timer.print_start()

    # Test the pipeline #
    print 'STEP 0: Testing the script and connection'
    test()
    print "-> Test OK !"
    timer.print_elapsed()

    # Do it #
    gis, length = all_gis()
    timer.print_elapsed()
    add_to_database(gis_to_records(gis))
    timer.print_elapsed()

    # End #
    print 'Done. Results are in "%s"' % os.path.abspath(current_dir)
    timer.print_end()
    timer.print_total_elapsed()