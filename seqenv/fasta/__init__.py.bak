# Built-in modules #
import gzip

# Internal modules #
from seqenv.common.cache import property_cached
from seqenv.common.autopaths import FilePath

# Third party modules #
import sh
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
class FASTA(FilePath):
    """A single FASTA file somewhere in the filesystem. You can read from it in
    several convenient ways. You can write to it in a automatically buffered way.
    There are several other things you can do with a FASTA file. Look at the class."""

    format    = 'fasta'
    extension = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()

    @property
    def gzipped(self): return True if self.path.endswith('gz') else False

    @property_cached
    def count(self):
        if self.gzipped: return int(sh.zgrep('-c', "^>", self.path, _ok_code=[0,1]))
        else: return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    def open(self, mode='r'):
        if self.gzipped: self.handle = gzip.open(self.path, mode)
        else: self.handle = open(self.path, mode)

    def close(self):
        if hasattr(self, 'buffer'):
            self.flush()
            del self.buffer
        self.handle.close()

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.format)

    def create(self):
        self.buffer = []
        self.buf_count = 0
        if not self.directory.exists: self.directory.create()
        self.open('w')
        return self

    def add_seq(self, seq):
        self.buffer.append(seq)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def add_str(self, seq, name=None, description=""):
        self.add_seq(SeqRecord(Seq(seq), id=name, description=description))

    def flush(self):
        for seq in self.buffer: SeqIO.write(seq, self.handle, self.format)
        self.buffer = []

    def write(self, reads):
        if not self.directory.exists: self.directory.create()
        self.open('w')
        SeqIO.write(reads, self.handle, self.format)
        self.close()

    @property_cached
    def ids(self):
        """All the sequences ids in a set"""
        return frozenset([seq.id for seq in self])

    def rename_sequences(self, new_fasta, mapping):
        """Given a new file path, will rename all sequences in the
        current fasta file using the mapping dictionary also provided."""
        assert isinstance(new_fasta, FASTA)
        new_fasta.create()
        for seq in self:
            new_name = mapping[seq.id]
            nucleotides = str(seq.seq)
            new_fasta.add_str(nucleotides, new_name)
        new_fasta.close()

    def extract_sequences(self, new_fasta, ids):
        """Will take all the sequences from the current file who's id appears in
        the ids given and place them in the new file path given."""
        assert isinstance(new_fasta, FASTA)
        new_fasta.create()
        for seq in self:
            if seq.id in ids: new_fasta.add_seq(seq)
        new_fasta.close()