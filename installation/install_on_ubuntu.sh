# To install (and test) seqenv starting from a clean Ubuntu 14.04.4 LTS

sudo apt-get install libgraphviz-dev libfreetype6-dev ncbi-blast+

cd ~
mkdir test
cd test

pip install --user virtualenvwrapper
source virtualenvwrapper_lazy.sh

export PYTHONPATH=""

mkvirtualenv kiwi
pip install numpy matplotlib
pip install seqenv
python -c "import seqenv; print seqenv.__version__; print seqenv"

cp ~/repos/seqenv/examples/minimal/test.fasta .
seqenv test.fasta --out_dir output

deactivate
rmvirtualenv kiwi