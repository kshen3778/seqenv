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
from seqenv.common.autopaths import FilePath
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
restore_path = FilePath("restore_gi_db.zip")
db_path      = FilePath("gi_db.sqlite3")

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
def restore_database(timer):
    print '-> Restoring database.'
    db_path.remove()
    restore_path.unzip_to(db_path)
    timer.print_elapsed()

###############################################################################
def pre_test(database, n=20):
    print '-'*50
    database.execute('SELECT * from data')
    tagger = Tagger()
    for i in range(n):
        gi, text, pubmed = database.cursor.next()
        matches = tagger.match(text)
        if not matches: continue
        envos = tuple(int(n[1][5:]) for m in matches for n in m[2])
        msg = 'GI: %s, SOURCE: "%s…", ENVOS: %s'
        msg = msg % (gi, text[:15], envos)
        print msg
    print '-'*50

###############################################################################
def run(database, timer, start_again=False):
    """Main part."""
    # STEP 0 #
    if start_again:
        print 'STEP 0: Starting over. Droping table.'
        database.execute('DROP TABLE IF EXISTS "isolation"')

    # STEP 1 #
    print 'STEP 1: Adding new table if not exists.'
    query = """
    CREATE TABLE IF NOT EXISTS "isolation"
    (
      "id"     INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
      "source" TEXT NOT NULL,
      "envos"  BLOB NOT NULL
    );"""
    database.execute(query)
    timer.print_elapsed()

    # STEP 2 #
    print 'STEP 2: Loading all distinct isolation sources in RAM (~1min).'
    def gen_sources():
        command = 'SELECT DISTINCT "source" FROM "data" order by min(rowid);'
        database.execute(command)
        for x in database.cursor: yield x[0]
    all_sources = list(gen_sources())
    timer.print_elapsed()
    print "Total sources: %i" % len(all_sources)

    # STEP 3A #
    print 'STEP 3A: Adding all isolation sources and tags in the new table (~2min).'
    tagger = Tagger()
    def gen_rows(sources):
        for i in trange(len(sources)):
            text    = sources[i]
            matches = tagger.match(text)
            if not matches: continue
            envos   = (int(n[1][5:]) for m in matches for n in m[2])
            blob    = marshal.dumps(tuple(envos))
            yield text, blob
    database.add(gen_rows(all_sources), 'isolation', ('source', 'envos'))
    timer.print_elapsed()
    total = database.count_entries('isolation')
    print "Total sources that had at least one match: %i" % total

    # STEP 3B #
    print 'STEP 3B: Indexing.'
    database.index('isolation', 'source')
    timer.print_elapsed()

    # STEP 4 #
    print 'STEP 4: Adding the key column for isolation source texts.'
    for gi, text, pubmed in database:
        result = database.get_entry(text, 'source', 'isolation')
        if not result: database.
        num, source, envos = result
        envos = marshal.loads(envos)
        print num, source, envos
        break
    timer.print_elapsed()

    # STEP 5 #
    print 'STEP 5: Removing the source text column and removing matchless entries.'
    for gi, text, pubmed in database:
        result = database.get_entry(text, 'source', 'isolation')
        if not result: database.
        num, source, envos = result
        envos = marshal.loads(envos)
        print num, source, envos
        break
    timer.print_elapsed()

    # STEP 6 #
    print 'STEP 6: Vacuuming.'
    database.vacuum()
    timer.print_elapsed()

    # End messages #
    print 'Done. Results are in "%s"' % database.path
    timer.print_end()
    timer.print_total_elapsed()

###############################################################################
def post_test(database, n=20):
    print '-'*50
    database.execute('SELECT * from "isolation"')
    for i in range(n/2):
        num, source, envos = database.cursor.next()
        envos = marshal.loads(envos)
        msg = 'NUM: %s, SOURCE: "%s…", ENVOS: %s'
        msg = msg % (num, source[:15], envos)
        print msg
    print '-'*50
    pass
    print '-'*50

###############################################################################
if __name__ == '__main__':
    print "*** Adding tagger results to database (pid %i) ***" % os.getpid()
    # Time it #
    timer = Timer()
    timer.print_start()
    # Unzip again #
    restore_database(timer)
    # Do it #
    with Database(db_path, text_fact=bytes) as db:
        pre_test(db)
        run(db, timer)
        post_test(db)