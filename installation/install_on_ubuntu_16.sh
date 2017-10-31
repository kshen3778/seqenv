# To install (and test) seqenv starting from a clean Ubuntu 16.04.3 LTS
# In this case the 64-bit version running inside VirtualBox

sudo apt-get install git pip libgraphviz-dev libfreetype6-dev ncbi-blast+

pip install --upgrade pip
pip install numpy matplotlib
pip install --user seqenv

mkdir ~/repos/
cd ~/repos/
git clone https://github.com/xapple/seqenv.git
cd ~/repos/seqenv/examples/minimal/

wget -O - ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz | gunzip -c | head -n 5000000 > fake_nt.fasta
makeblastdb -in fake_nt.fasta -dbtype nucl

seqenv test.fasta --out_dir output --search_db fake_nt.fasta