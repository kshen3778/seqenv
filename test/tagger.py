"""
===============
Test the tagger
===============
"""

# Modules #
import tagger, os

# Script #
t = tagger.Tagger()
home = os.environ['HOME'] + '/'
data_dir = home + "repos/seqenv/data/"
t.LoadNames(data_dir + 'envo_entities.tsv', data_dir + 'envo_names.tsv')
t.LoadGlobal(data_dir + 'envo_global.tsv')
t.GetMatches("oceanic soil", "", [-27])