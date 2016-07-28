# To install (and test) seqenv on macOS 10.11

brew install blast
brew install graphviz --with-bindings

cd ~
mkdir test
cd test

pip install --user virtualenvwrapper
source virtualenvwrapper_lazy.sh

export PYTHONPATH=""

mkvirtualenv kiwi
pip install numpy
pip install --global-option=build_ext --global-option="-I/usr/local/Cellar/graphviz/2.38.0/include/"  --global-option="-L/usr/local/Cellar/graphviz/2.38.0/lib/" pygraphviz

pip install seqenv
python -c "import seqenv; print seqenv.__version__; print seqenv"

wget -O - ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz | gunzip -c | head -n 10000 > fake_nt.fasta
makeblastdb -in fake_nt.fasta -dbtype nucl
printf ">test\nAAATAAGGACTTGTATGAATGGCCACACGAGGGTTTTACTGTCTCTTACTTTTAATCAGTGAAATTGACCTCCCCGTGAAGAGGCGGGGATAATAAAATAAGACGAGAAGACCCTATGGAGCTTTAATTAATCAACTCAAAAATCAGAAAACAATACCACTAAGGGATAACAGA" > seq.fasta
seqenv seq.fasta --out_dir output --search_db fake_nt.fasta

deactivate
rmvirtualenv kiwi