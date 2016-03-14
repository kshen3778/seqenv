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
import os, inspect, marshal

# Internal modules #
from seqenv.common.timer import Timer
from seqenv.common.database import Database
import seqenv.tagger as api

# Third party modules #
from tqdm import tqdm

# Get the directory of this script #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename))  + '/'
default_data_dir = current_dir + '../data_envo/'

# The default location of the database #
db_path = "gi_db.sqlite3"

###############################################################################
class Tagger(object):
    """Interface to the C-coded tagger. Randomly segfaults last
    time I tried it on a new machine."""

    def __init__(self, entities=None, names=None, globs=None, data_dir=None):
        # Defaults #
        if data_dir is None: data_dir = default_data_dir
        if entities is None: entities = data_dir + 'envo_entities.tsv'
        if names is None:    names    = data_dir + 'envo_names.tsv'
        if globs is None:    globs    = data_dir + 'envo_global.tsv'
        # Make an instance of the API #
        self.api = api.Tagger()
        self.api.LoadNames(entities, names)
        self.api.LoadGlobal(globs)

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

    # STEP 1 #
    print 'STEP 1: Adding new table if not exists.'
    database = Database(db_path)
    query = """
    CREATE TABLE IF NOT EXISTS isolation
    (
      source INTEGER AUTO_INCREMENT PRIMARY KEY,
      text TEXT,
      envo BLOB
    );
    """
    database.execute(query)
    timer.print_elapsed()

    # STEP 2 #
    print 'STEP 2: Loading all distinct isolation sources in RAM.'
    def gen_sources():
        query = "SELECT DISTINCT source FROM data;"
        for x in tqdm(database.execute(query)): yield x
    all_sources = list(gen_sources())
    timer.print_elapsed()

    # STEP 3 #
    print 'STEP 3: Finding resume point.'
    if database.count_entries('isolation') != 0:
        last = database.get_last('isolation')
        index = all_sources.index(last)
    else:
        index = 0
    timer.print_elapsed()

    # STEP 4 #
    print 'STEP 4: Adding all isolation sources and tags in the new table.'
    tagger = Tagger()
    def gen_rows(sources, i):
        for i in tqdm(xrange(index,len(sources))):
            text = sources[i]
            matches = tagger.match(text)
            blob = marshal.dumps((m for m in matches))
            yield i, text, blob
    database.add(gen_rows(all_sources, index))
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