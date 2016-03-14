#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
A script to add the results from the tagger to the database in a static fashion.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ ./add_tagger_results.py
"""

# Built-in modules #
import os, inspect

# Internal modules #
from seqenv import Analysis, repos_dir
from seqenv.common.timer import Timer
from seqenv.common.database import Database
import seqenv.tagger as api

# Third party modules #

# Get the directory of this script #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename))  + '/'
default_data_dir = current_dir + '../data_envo/'

# The default location of the database #
db_path = "gi_db.sqlite3"

###############################################################################
def modify_the_database():
    """Run this only once."""
    # Start timer #
    timer = Timer()
    timer.print_start()
    # Do it #
    print 'STEP 1: Adding new table.'
    database = Database(db_path)
    query = """
    CREATE TABLE isolation
    (
      source INTEGER AUTO_INCREMENT PRIMARY KEY,
      text TEXT
      tagged TEXT
    );
    """
    database.execute(query)
    timer.print_elapsed()
    # Do it #
    print 'STEP 2: Adding all isolation sources.'
    analysis = Analysis(repos_dir + "examples/minimal/test.fasta")
    analysis.add
    timer.print_elapsed()
    # End messages #
    timer.print_elapsed()
    print 'Done. Results are in "%s"' % 'a'
    timer.print_end()
    timer.print_total_elapsed()

###############################################################################
class Tagger(object):
    """Interface to the C-coded tagger. Randomly segfaults last
    time I tried it on a new machine."""

    def __init__(self, entities=None, names=None, globals=None, data_dir=None):
        # Defaults #
        if data_dir is None: data_dir = default_data_dir
        if entities is None: entities = data_dir + 'envo_entities.tsv'
        if names is None:    names = data_dir + 'envo_names.tsv'
        if globals is None:  globals = data_dir + 'envo_global.tsv'
        # Make an instance of the API #
        self.api = api.Tagger()
        self.api.LoadNames(entities, names)
        self.api.LoadGlobal(globals)

    def match(self, text):
        """When you call `t.GetMatches(text, "", [-27])` you get a list back.
        The second argument can be left empty in our case (per document blacklisting)
        The result is something like:
        - [(53, 62, ((-27, 'ENVO:00002001'),)), (64, 69, ((-27, 'ENVO:00002044'),))]
        The number -27 is ENVO terms, -26 could be tissues, etc.
        Sometimes, one word can link to two concepts such as with the word 'marine'."""
        return self.api.GetMatches(text, "", [-27])

###############################################################################
def run():
    """Run this script."""
    # Start timer #
    timer = Timer()
    timer.print_start()
    # Do it #
    print 'STEP 1: Finding resume point.'
    database = Database(db_path)
    database.cursor("")
    timer.print_elapsed()
    # Do it #
    print 'STEP 2: Adding tagger results.'
    timer.print_elapsed()
    # End messages #
    timer.print_elapsed()
    print 'Done. Results are in "%s"' % 'a'
    timer.print_end()
    timer.print_total_elapsed()

###############################################################################
if __name__ == '__main__':
    print "*** Adding tagger results to database (pid %i) ***" % os.getpid()
    pass