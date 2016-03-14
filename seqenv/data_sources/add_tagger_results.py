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
import os, inspect, marshal, unicodedata

# Internal modules #
from seqenv.common.timer import Timer
from seqenv.common.database import Database
import tagger as tagger_api

# Third party modules #
from tqdm import trange

# Get the directory of this script #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename))  + '/'
default_data_dir = current_dir + '../data_envo/'

# The default location of the database #
db_path = "gi_db.sqlite3"

###############################################################################
class Tagger(object):
    """Interface to the C-coded tagger. Randomly segfaults last
    time I recomplied it and then tried it on a new machine."""

    def __init__(self, entities=None, names=None, globs=None, data_dir=None):
        # Defaults #
        if data_dir is None: data_dir = default_data_dir
        if entities is None: entities = data_dir + 'envo_entities.tsv'
        if names is None:    names    = data_dir + 'envo_names.tsv'
        if globs is None:    globs    = data_dir + 'envo_global.tsv'
        # Make an instance of the API #
        self.api = tagger_api.Tagger()
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
def run(database, start_again=True):
    """Run this script."""
    # Start timer #
    timer = Timer()
    timer.print_start()

    # STEP 0 #
    if start_again: database.execute('DROP TABLE IF EXISTS "isolation"')

    # STEP 1 #
    print 'STEP 1: Adding new table if not exists.'
    query = """
    CREATE TABLE IF NOT EXISTS "isolation"
    (
      "source" INTEGER AUTO_INCREMENT PRIMARY KEY,
      "text" TEXT,
      "envos" BLOB
    );"""
    database.execute(query)
    timer.print_elapsed()

    # STEP 2 #
    print 'STEP 2: Loading all distinct isolation sources in RAM (1min).'
    def gen_sources():
        command = 'SELECT DISTINCT "source" FROM "data";'
        database.execute(command)
        for x in tqdm(database.cursor): yield x[0]
    all_sources = list(gen_sources())
    timer.print_elapsed()
    print "Total sources: %i" % len(all_sources)

    # STEP 3 #
    print 'STEP 3: Finding resume point.'
    if database.count_entries('isolation') != 0:
        last = database.get_last('isolation')
        index = all_sources.index(last) + 1
    else:
        index = 0
    timer.print_elapsed()

    # STEP 4 #
    print 'STEP 4: Adding all isolation sources and tags in the new table.'
    tagger = Tagger()
    def gen_rows(sources, i):
        for i in trange(index,len(sources)):
            text    = sources[i]
            ascii   = unicodedata.normalize('NFKD', text).encode('ascii','ignore')
            matches = tagger.match(ascii)
            if not matches: continue
            envos   = (int(n[1][5:]) for m in matches for n in m[2])
            blob    = marshal.dumps(tuple(envos))
            yield text, blob
    database.add(gen_rows(all_sources, index), 'isolation', ('text', 'envos'))
    timer.print_elapsed()
    total = database.count_entries('isolation')
    print "Total sources that had at least one match: %i" % total

    # End messages #
    print 'Done. Results are in "%s"' % database.path
    timer.print_end()
    timer.print_total_elapsed()

###############################################################################
if __name__ == '__main__':
    print "*** Adding tagger results to database (pid %i) ***" % os.getpid()
    with Database(db_path) as db: run(db)