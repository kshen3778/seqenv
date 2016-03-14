"""
Test the tagger
"""

# Modules #
import tagger, os, inspect

# Get directory #
filename = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'
data_dir = current_dir +  "../seqenv/data_envo/"

# Script #
t = tagger.Tagger()
t.LoadNames(data_dir + 'envo_entities.tsv', data_dir + 'envo_names.tsv')
t.LoadGlobal(data_dir + 'envo_global.tsv')
results = t.GetMatches("oceanic soil", "", [-27])
print results