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
from tqdm import trange, tqdm

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
def make_test_database(database):
    command = """
    DELETE FROM "data"
        WHERE id NOT IN (SELECT id
                         FROM "data"
                         LIMIT 10000);"""
    database.execute(command)

###############################################################################
def restore_database(timer):
    print '-> Restoring database (~5sec).'
    db_path.remove()
    restore_path.unzip_to(db_path)
    timer.print_elapsed()

###############################################################################
def pre_test(database, n=20):
    print '-'*25 + "PRE TEST" + '-'*25
    database.execute('SELECT * from data')
    tagger = Tagger()
    ids = []
    for i in range(n):
        gi, text, pubmed = database.cursor.next()
        matches = tagger.match(text)
        if not matches: continue
        envos = tuple(int(n[1][5:]) for m in matches for n in m[2])
        msg = 'GI: %s, SOURCE: "%s…", ENVOS: %s'
        msg = msg % (gi, text[:15], envos)
        ids.append(gi)
        print msg
    print '-'*50
    return ids

###############################################################################
def run(database, timer):
    """Main part."""
    # STEP 1 #
    print 'STEP 1: Adding the isolation table.'
    query = """
    CREATE TABLE "isolation"
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

    # STEP 3 #
    print 'STEP 3: Adding all isolation sources and tags in the new table (~1min).'
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

    # STEP 4 #
    print 'STEP 4: Indexing on source text in the isolation table.'
    database.index(table='isolation', column='source')
    timer.print_elapsed()

    # STEP 5 #
    print 'STEP 5: Uniquifying Gids in the old table. (~1min)'
    query = 'DELETE from "data" WHERE rowid not in (select min(rowid) from data group by id);'
    database.execute(query)
    timer.print_elapsed()

    # STEP 6 #
    query = """
    CREATE TABLE "gi"
    (
      "id"      INTEGER PRIMARY KEY NOT NULL,
      "isokey"  INTEGER NOT NULL REFERENCES "isolation"("id"),
      "pubmed"  INTEGER
    );"""
    database.execute(query)
    print 'STEP 6A: Putting the new table in RAM. (~3min)'
    def gen_rows(database):
        for gi, text, pubmed in tqdm(database):
            isolation = database.get_entry(text, 'source', 'isolation')
            if isolation: yield gi, isolation[0], pubmed
    rows = tuple(gen_rows(database))
    print 'STEP 6B: Filling the new table into the database. (~30min)'
    database.add(tqdm(rows), 'gi')
    timer.print_elapsed()

    # STEP 5 #
    print 'STEP 5: Delete the old table, index and make a view.'
    command = 'DROP INDEX "main_index";'
    database.execute(command)
    command = 'DROP TABLE "data";'
    database.execute(command)
    database.index(table='gi', column='id')
    command   = 'CREATE VIEW "data" AS'
    command  += ' SELECT data.id, %s from "gi";'
    subselect = '(SELECT source from isolation where gi.isokey=isolation.id)'
    database.execute(command % subselect)
    timer.print_elapsed()

    # STEP 6 #
    print 'STEP 6: Rebuild the entire database (vacuuming). (~xmin)'
    database.vacuum()
    timer.print_elapsed()

    # End messages #
    print 'Done. Results are in "%s"' % database.path
    timer.print_end()
    timer.print_total_elapsed()

###############################################################################
def post_test(database, ids, n=20):
    print '-'*25 + "POST TEST 1" + '-'*25
    database.execute('SELECT rowid, * from "isolation"')
    for i in range(n/2):
        rowid, key, source, envos = database.cursor.next()
        envos = marshal.loads(envos)
        msg = 'ID: %s, KEY: %s, SOURCE: "%s…", ENVOS: %s'
        msg = msg % (rowid, key, source[:15], envos)
        print msg
    print '-'*25 + "POST TEST 2" + '-'*25
    command = """
    SELECT rowid, * from "data"
    WHERE rowid in (%s)
    ORDER BY CASE rowid
    %s
    END
    """
    ordered = ','.join(map(str,ids))
    rowids  = '\n'.join("WHEN '%s' THEN %s" % (row,i) for i,row in enumerate(ids))
    command = command % (ordered, rowids)
    # This could have worked but sqlite3 was too old on the server
    # ORDER BY instr(',%s,', ',' || id || ',')
    database.execute(command)
    for i in range(len(ids)):
        rowid, gi, isokey, pubmed = database.cursor.next()
        key, source, envos = database.get_entry(isokey, 'id', 'isolation')
        envos = marshal.loads(envos)
        msg = 'ID: %s, GI: %s, SOURCE: "%s…", ENVOS: %s'
        msg = msg % (rowid, gi, source[:15], envos)
        print msg
    print '-'*70

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
        ids = pre_test(db)
        run(db, timer)
        post_test(db, ids)
    # Size #
    print "Final size is %s" % db_path.size