#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
A script to parse and download information from the NCBI NT database
useful for the seqenv project. We want isolation sources and pubmed_ids
for all GI numbers found in the NT database.

See https://biosupport.se/p/719/
And https://www.biostars.org/p/153833/

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ ./generate.py
"""

# Built-in modules #
import os, time, inspect, urllib2, datetime
from itertools import islice
from collections import OrderedDict, defaultdict
from socket import error as SocketError

# Internal modules #
from seqenv.common           import GenWithLength, pretty_now
from seqenv.common.timer     import Timer
from seqenv.common.autopaths import FilePath
from seqenv.common.database  import Database
from seqenv.common.tmpstuff  import new_temp_file

# Third party modules #
from Bio import Entrez
from tqdm import tqdm
from shell_command import shell_output
from Bio.Entrez.Parser import CorruptedXMLError, ValidationError

# Get the directory of this script #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'

###############################################################################
class GenerateGIs(FilePath):
    """Takes care of generating a suite of GI numbers."""

    length_cutoff = 4000

    def retrieve_from_nt(self):
        """Get all GI numbers with their length from local NT database.
        Then filter out the ones that are too long."""
        temp_file = new_temp_file()
        shell_output("blastdbcmd -db nt -entry all -outfmt '%g %l' > " + temp_file)
        with_len = (map(int,line.strip('\n').split()) for line in temp_file)
        self.writelines(str(gid) + '\n' for gid, length in with_len if length < self.length_cutoff)
        temp_file.remove()

    def check_retrieved(self):
        """Check that we have retrieved the GI numbers, else do it."""
        if self.count_bytes < 1: self.retrieve_from_nt()
        else: print '-> All GI numbers already found in file "%s", skipping STEP 1.' % self.filename
        print '-> Got %i GI numbers.' % len(self)

    def all(self):
        """Iterate over all GI numbers."""
        return GenWithLength(self.lines_int, len(self))

    @property
    def lines_int(self):
        """An iterator on the lines of the file as integers"""
        for x in self: yield int(x.rstrip('\n'))

    def start_at(self, last_gi):
        """Iterate over all GI numbers starting at a given GI number."""
        generator = self.lines_int
        n = 0
        while True:
            gi = generator.next()
            n += 1
            if gi == last_gi: break
        else:
            raise Exception("Didn't find the GI '%s' in the all_gis.txt file" % last_gi)
        return GenWithLength(generator, len(self) - n)

def test_generate_gis():
    """A small generator that yields a couple test GI numbers.
    Useful for test cases."""
    gis = ['6451693', '127', '76365841', '22506766', '389043336',
           '497', '429143984', '264670502', '74268401', '324498487']
    return GenWithLength((gi for gi in gis), len(gis))

###############################################################################
class QueryNCBI(object):
    """Takes care of querying NCBI and restarting when errors occur"""

    email = "I don't know who will be running this script"

    def __init__(self, at_a_time=100, delay=0.2, progress=True):
        """Give the number of records to retrieve at a time as well
        as the delay to wait between requests"""
        # Base parameters #
        self.at_a_time = int(at_a_time)
        self.delay     = float(delay)
        # Display progress #
        if progress is True:    self.progress = tqdm
        elif progress is False: self.progress = lambda x:x
        else:                   self.progress = progress
        # Set the email #
        Entrez.email = self.email
        # Error logging #
        self.logger = Logger(current_dir + 'error_log.txt')
        # Other #
        self.chunk = []

    def get(self, gi_nums):
        """Download information from NCBI in batch mode.
        Return `isolation_source` and `pubmed_id for an entry if it has an
        isolation_source"""
        for i in self.progress(xrange(0, len(gi_nums), self.at_a_time)):
            self.chunk = tuple(islice(gi_nums, 0, self.at_a_time))
            records    = self.chunk_to_records(self.chunk)
            sources    = map(self.record_to_source, records)
            pubmed_ids = map(self.record_to_pubmed_id, records)
            has_source = tuple(i for i in xrange(len(self.chunk)) if sources[i])
            # Generate the result #
            if not has_source: continue
            yield ((self.chunk[i], sources[i], pubmed_ids[i]) for i in has_source)

    def chunk_to_records(self, chunk):
        """Download from NCBI until it works. Will restart until reaching the python
        recursion limit. We don't want to get our IP banned from NCBI so we have
        a little pause at every function call."""
        time.sleep(self.delay)
        # Check that none of the sequences are above a certain threshold in length ? #
        pass
        # The web connection #
        try:
            response = Entrez.efetch(db="nucleotide", id=map(str,chunk), retmode="xml")
        except urllib2.HTTPError:
            self.logger.report_error('http', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        except urllib2.URLError:
            self.logger.report_error('url', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        # The parsing xml #
        try:
            result = list(Entrez.parse(response))
        except CorruptedXMLError:
            self.logger.report_error('xml', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        except ValidationError:
            self.logger.report_error('parse', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        except SocketError:
            self.logger.report_error('socket', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        # Sometimes it's one short xD #
        if len(result) != len(chunk):
            self.logger.report_error('length', chunk)
            time.sleep(5)
            return self.chunk_to_records(chunk)
        # Return the result #
        assert len(result) == len(chunk)
        return result

    def record_to_source(self, record):
        qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
        for qualifier in qualifiers:
            if qualifier['GBQualifier_name'] == 'isolation_source':
                return qualifier['GBQualifier_value']

    def record_to_pubmed_id(self, record):
        if 'GBSeq_references' not in record: return None
        references = record['GBSeq_references'][0]
        if 'GBReference_pubmed' not in references: return None
        return int(references.get('GBReference_pubmed'))

    def get_all_lengths(self, chunk):
        """Given a list of GIDs return the length of each sequence"""
        gid_to_len = lambda x: Entrez.parse(Entrez.esummary(db="nucleotide", id=x)).next()['Length']
        return map(gid_to_len, chunk)

###############################################################################
class Logger(FilePath):
    """Takes care of error logging for the NCBI worker."""

    def __init__(self, path):
        self.path = self.clean_path(path)
        self.errors = defaultdict(list)

    def report_error(self, kind, *args, **kwargs):
        self.errors[kind] += [datetime.datetime.now()]
        self.write(self.message)

    @property
    def message(self):
        message = ''
        message += "Last updated:       %s\n" % pretty_now()
        message += "HTTP errors:        %s\n" % len(self.errors['http'])
        message += "URL errors:         %s\n" % len(self.errors['url'])
        message += "XML errors:         %s\n" % len(self.errors['xml'])
        message += "Parse errors:       %s\n" % len(self.errors['parse'])
        message += "Socket errors:      %s\n" % len(self.errors['socket'])
        message += "Freak length error: %s\n" % len(self.errors['length'])
        return message

    def print_error_log(self): print self.message

###############################################################################
# Objects #
gi_generator = GenerateGIs(current_dir + 'small_gis.txt')
sqlite_db    = Database(current_dir + 'gi_db.sqlite3')
ncbi_worker  = QueryNCBI()

###############################################################################
def run():
    """Run this script."""
    # Start timer #
    timer = Timer()
    timer.print_start()
    # Get the numbers #
    print 'STEP 1: Get all GIs from the local nt database into a file (about 6h).'
    gi_generator.check_retrieved()
    timer.print_elapsed()
    # Do it #
    print 'STEP 2: Querying NCBI and writing to database (about 150h)'
    if sqlite_db.count_bytes < 1:
        print "-> The database doesn't exist. Creating it."
        keys = OrderedDict((('id','integer'), ('source','text'), ('pubmed','integer')))
        sqlite_db.create(keys)
        sqlite_db.add_by_steps(ncbi_worker.get(gi_generator.all()))
    else:
        print '-> The database seems to exist already. Attempting to restart where it was stopped.'
        last_gi = sqlite_db.last[0]
        print '-> Last GI number recorded was "%s".' % last_gi
        sqlite_db.add_by_steps(ncbi_worker.get(gi_generator.start_at(last_gi)))
    # Make the index #
    print 'STEP 3: Creating database index.'
    sqlite_db.index()
    # End messages #
    timer.print_elapsed()
    print 'Done. Results are in "%s"' % sqlite_db
    timer.print_end()
    timer.print_total_elapsed()
    ncbi_worker.logger.print_error_log()

###############################################################################
def test():
    """Test this script."""
    # Verbose #
    print 'STEP 0: Testing the script and connection'
    # The result we should get #
    template_result = """
    76365841,Everglades wetlands,16907754
    389043336,lake water at 5 m depth during dry season,None
    429143984,downstream along river bank,None
    264670502,aphotic layer; anoxic zone; tucurui hydroeletric power plant reservoir,None
    324498487,bacterioplankton sample from lake,None"""
    # Make it pretty #
    template_result = '\n'.join(l.lstrip(' ') for l in template_result.split('\n') if l)
    # The result we got #
    test_ncbi_worker = QueryNCBI(at_a_time=10, progress=False)
    results = test_ncbi_worker.get(test_generate_gis()).next()
    # Function #
    def result_to_lines(results):
        for gi, source, pubmed in results:
            yield str(gi) + ',' + source + ',' + str(pubmed) + '\n'
    # Check #
    text_version = ''.join(result_to_lines(results))
    assert text_version.strip('\n') == template_result.strip('\n')
    print "-> Test OK !"

###############################################################################
def validate():
    """Check that the database is good."""
    pass

###############################################################################
if __name__ == '__main__':
    print "*** Generating the GI database (pid %i) ***" % os.getpid()
    test()
    run()
    validate()