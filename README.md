[![changelog](http://allmychanges.com/p/python/seqenv/badge/)](http://allmychanges.com/p/python/seqenv/?utm_source=badge) [![PyPI version](https://badge.fury.io/py/seqenv.svg)](https://badge.fury.io/py/seqenv)

# `seqenv` version 1.3.0
* Assign environment ontology (EnvO) terms to short DNA sequences.
* All code written by [Lucas Sinclair](http://envonautics.com/#lucas).
* Publication at: https://peerj.com/articles/2690/

### Dependencies
* You need to have a copy of the NCBI nucleotide data base (called `nt`) installed locally as well as the `blastn` executable in your `$PATH`. So that BLAST finds the NT database, you can edit your `~/.ncbirc` file. The `blastn` executable should be in your `$PATH`.
* You need to have `graphviz` installed and the C libraries should be accessible. This is usually resolved by typing `apt-get install libgraphviz-dev`. The `dot` executable should be in your `$PATH`.
* The project also depends on some other python modules such as `biopython` and `pandas`. Happily, these will be installed automatically when calling the `pip` command below.

### Installing
To install `seqenv` onto your machine, use the python package manager:

    $ pip install seqenv

You might be installing this onto a computer server which you don't own and thus don't have sufficient privileges. In that case you can install everything in your home directory like this:

    $ pip install --user seqenv

If this still doesn't work, you might be missing the `pip` program on your system or the correct version of Python (any version `2.7.x`). You can get both of these things by using using this little project: https://github.com/yyuu/pyenv

### Usage
Once that is done, you can start processing FASTA files from the command line. For using the default parameters you can just type:

    $ seqenv sequences.fasta

We will then assume that you have inputed 16S sequences. To modify the database or input a different type of sequences:

    $ seqenv sequences.fasta --seqtype prot --search_db nr

To modify the minimum identity in the similarity search, use the following:

    $ seqenv sequences.fasta --min_identity 0.97

If you have abundance data you would like to add to your analysis you can specify it like this in a TSV file:

    $ seqenv sequences.fasta --abundances counts.tsv

### All parameters
Several other options are possible. Here is a list describing them all:

   * `--seq_type`: Sequence type `nucl` or `prot`, for nucleotides or amino acids respectively (Default: `nucl`).
   * `--search_algo`: Search algorithm. Either `blast` or `vsearch` (Default: `blast`).
   * `--search_db`: The database to search against (Default: `nt`). You can specify the full path or make a `~/.ncbirc` file.
   * `--normalization`: Can be either of `flat`, `ui` or `upui`. This option defaults to `flat`.
     * If you choose `flat`, we will count every isolation source independently,
       even if the same text appears several times for the same input sequence.
     * If you choose `ui`, standing for unique isolation, we will count every
       identical isolation source only once within the same input sequence.
     * If you choose `upui`, standing for unique isolation and unique pubmed-ID,
       we will uniquify the counts based on the text entry of the isolation
       sources as well as on the pubmed identifiers from which the GI obtained.
   * `--proportional`: Should we divide the counts of every input sequence by the number of envo terms that were associated to it. Defaults to `True`.
   * `--backtracking`: For every term identified by the tagger, we will propagate frequency counts up the acyclic directed graph described by the ontology. Defaults to `False`.
   * `--restrict`: Restrict the output to the descendants of just one ENVO term. This removes all other terms that are not reachable through the given node. For instance you could specify: `ENVO:00010483` (Disabled by default)
   * `--num_threads`: Number of cores to use (Defaults to the total number of cores). Use `1` for non-parallel processing.
   * `--out_dir`: The output directory in which to store the result and intermediary files. Defaults to the same directory as the input file.
   * `--min_identity`: Minimum identity in similarity search (Default: `0.97`). Note: not available when using `blastp`.
   * `--e_value`: Minimum e-value in similarity search (Default: `0.0001`).
   * `--max_targets`: Maximum number of reference matches in the similarity search (Default: `10`).
   * `--min_coverage`: Minimum query coverage in similarity search (Default: `0.97`).
   * `--abundances`: Abundances file as TSV with OTUs as rows and sample names as columns (Default: None).
   * `--N`: If abundances are given, pick only the top N sequences (Disabled by default).

### Why make this ?
The continuous drop in the associated costs combined with the increased efficiency of the latest high-throughput sequencing technologies has resulted in an unprecedented growth in sequencing projects. Ongoing endeavors such as the [Earth Microbiome Project](http://www.earthmicrobiome.org) and the [Ocean Sampling Day](http://www.microb3.eu/osd) are transcending national boundaries and are attempting to characterize the global microbial taxonomic and functional diversity for the benefit of mankind. The collection of sequencing information generated by such efforts is vital to shed light on the ecological features and the processes characterizing different ecosystems, yet, the full knowledge discovery potential can only be unleashed if the associated meta data is also exploited to extract hidden patterns. For example, with the majority of genomes submitted to NCBI, there is an associated PubMed publication and in some cases there is a GenBank field called "isolation sources" that contains rich environmental information.

With the advances in community-generated standards and the adherence to recommended annotation guidelines such as those of [MIxS](http://gensc.org/gc_wiki/index.php/MIxS) of the Genomics Standards Consortium, it is now feasible to support intelligent queries and automated inference on such text resources.

The [Environmental Ontology](http://environmentontology.org/) (or EnvO) will be a critical part of this approach as it gives the ontology for the concise, controlled description of environments. It thus provides structured and controlled vocabulary for the unified meta data annotation, and also serves as a source for naming environmental information. Thus, we have developed the `seqenv` pipeline capable of annotating sequences with environment descriptive terms occurring within their records and/or in relevant literature.

The `seqenv` pipeline can be applied to any set of nucleotide or protein sequences. Annotation of metagenomic samples, in particular 16S rRNA sequences is also supported.

The pipeline has already been applied to a range of datasets (e.g Greek lagoon, Swedish lake/river, African and Asian pitlatrine datasets, Black Sea sediment sample datasets have been processed).

### What does it do exactly ?
 Given a set of DNA sequences, `seqenv` first retrieves highly similar sequences from public repositories (e.g. NCBI GenBank) using BLAST or similar algorithm. Subsequently, from each of these homologous records, text fields carrying environmental context information such as the reference title and the **isolation source** field found in the metadata are extracted. Once the relevant pieces of text from each matching sequence have been gathered, they are processed by a text mining module capable of identifying any EnvO terms they contain (e.g. the word "glacier", or "pelagic", "forest", etc.). The identified EnvO terms along with their frequencies of occurrence are then subjected to multivariate statistics, producing matrices relating your samples to their putative sources as well as other useful outputs.

### Pipeline overview
The publication contains more information of course, but here is a schematic overview of what happens inside `seqenv`:

[![seqenv](documentation/frequencies.png)](documentation/frequencies.png)

### Tutorial
We will first run `seqenv` on a 16S rRNA dataset using ***isolation sources*** as a text source. Here, `abundance.tsv` is a species abundance file (97% OTUs) processed through [`illumitag`](https://github.com/limno/illumitag) software and `centers.fasta` contains the corresponding sequences for the OTUs.

~~~
$ ls
abundance.tsv
centers.fasta

$ seqenv centers.fasta --abundances abundance.tsv --seq_type nucl --out_dir output --N 1000 --min_identity 0.99
~~~

The output you will receive should look something like this:

~~~
seqenv version 1.2.0 (pid 52169)
Start at: 2016-03-02 00:22:09.727377
--> STEP 1: Parse the input FASTA file.
Elapsed time: 0:00:00.005811
Using: output/renamed.fasta
--> STEP 2: Similarity search against the 'nt' database with 5 processes
Elapsed time: 0:02:11.215829
--> STEP 3: Filter out bad hits from the search results
Elapsed time: 0:00:00.002071
--> STEP 4: Parsing the search results
Elapsed time: 0:00:00.002099
--> STEP 5: Setting up the SQLite3 database connection.
Elapsed time: 0:00:00.054077
Got 81 GI hits and 65 of them had one for more EnvO terms associated.
--> STEP 6: Computing EnvO term frequencies.
Elapsed time: 0:00:00.721455
------------
Success. Outputs are in 'output/'
End at: 2016-03-02 00:24:22.504485
Total elapsed time: 0:02:12.777297
~~~

Once the pipeline has finished processing, you will have the following contents in the output folder:

~~~
$ ls output/
list_concepts_found.tsv  samples_to_names.tsv  seq_to_names.tsv   top_seqs.fasta.parts
renamed.fasta            seq_to_concepts.tsv   top_seqs.blastout
samples.biom             seq_to_gis.pickle     top_seqs.fasta
~~~

The most interesting files are probably:

* `list_concepts_found.tsv` links every OTU to all its relevant BLAST hits and linked ENVO terms.
* `seq_to_names.tsv` a matrix linking every OTU to its "composition" in terms of ENVO identifiers translated to readable names.
* `samples_to_names.tsv` if an abundance file was provided, this is a a matrix linking every one of your samples to its "composition" in terms of ENVO identifiers translated to readable names.
* `graphviz/` directory containing ontology graphs for everyone of the inputed sequences, such as in the following example:

[![seqenv](documentation/ontology_graph.png)](documentation/ontology_graph.png)

### Acknowledgments
`seqenv` was originally conceived in the following hackathons supported by European Union's Earth System Science and Environmental Management ES1103 COST Action ("[Microbial ecology & the earth system: collaborating for insight and success with the new generation of sequencing tools](http://www.cost.eu/domains_actions/essem/Actions/ES1103)"):

- **From Signals to Environmentally Tagged Sequences** (Ref: ECOST-MEETING-ES1103-050912-018418), September 27th-29th 2012, Hellenic Center for Marine Research, Crete, Greece.
- **From Signals to Environmentally Tagged Sequences II** (Ref: ECOST-MEETING-ES1103-100613-031037), June 10th-13th 2013, Hellenic Center for Marine Research, Crete, Greece.
- **From Signals to Environmentally Tagged Sequences III** (Ref: ECOST-MEETING-ES1103-220914-047036), September 22nd-25th 2014, Hellenic Center for Marine Research, Crete, Greece.

This work would not have been possible without the advice and support of many people who attended the hackathons, in alphabetical order:

- Simon Berger (simon.berger@h-its.org)
- Alica Chroňáková (alicach@upb.cas.cz)
- Lars Juhl Jensen (lars.juhl.jensen@cpr.ku.dk)
- Anastasis Oulas (oulas@hcmr.gr)
- [Evangelos Pafilis](http://epafilis.info/) (pafilis@hcmr.gr) [2]
- Christina Pavloudi (cpavloud@hcmr.gr)
- [Chris Quince](http://www2.warwick.ac.uk/fac/med/staff/cquince/) (c.quince@warwick.ac.uk) [3]
- Julia Schnetzer (jschnetz@mpi-bremen.de)
- [Lucas Sinclair](http://envonautics.com/#lucas) (lucas.sinclair@me.com) [1]
- Aaron Weimann (aaron.weimann@uni-duesseldorf.de)
- Ali Zeeshan Ijaz (alizeeshanijaz@gmail.com)
- [Umer Zeeshan Ijaz](http://userweb.eng.gla.ac.uk/umer.ijaz/) (umer.ijaz@glasgow.ac.uk) [3]

[1] Main developer
[2] Contact for correspondence
[3] Original idea

### News
* **August 2013**: Chris Quince presented a talk on `seqenv` at [STAMPS2013](https://stamps.mbl.edu/index.php/Main_Page). You can download the PDF of the presentation: [C Quince et. al., SeqEnv: Annotating sequences with environments (STAMPS 2013)](https://stamps.mbl.edu/images/4/44/Quince_SeqEnvSTAMPS2013.pdf)

### Publications that have used `seqenv`
* Bacterial diversity along a 2600 km river continuum. `doi:10.1111/1462-2920.12886`
* Can marine bacteria be recruited from freshwater sources and the air? `doi:10.1038/ismej.2014.89`

### Development installation
This chapter shows you how to get a development install that makes it easier to contribute to and to change `seqenv` to suit your needs. If you cannot get a functional installation set up, feel free to contact the authors.

##### Step 1: Cloning the repository
Here you will download a copy of the code and place it somewhere in your home directory.

    $ cd ~
    $ mkdir repos
    $ cd repos
    $ git clone https://github.com/xapple/seqenv.git

##### Step 2: Modify your search paths
Here you will edit your ``~/.bashrc`` or ``~/.bash_profile`` to add a reference to the code you just downloaded. You need to add these two lines:

    $ vim ~/.bash_profile
    export PYTHONPATH="$HOME/repos/seqenv/":$PYTHONPATH
    export PATH="$HOME/repos/seqenv/seqenv":$PATH

And finally source your profile file if you haven't already.

##### Step 3 (optional): Install your own version of python
Your system probably comes with a version of python installed. But the variations from system to system are too great to rely on any available python. We suggest to just install our own version in your home directory. Otherwise, make sure that you are using version 2.7.x of python.

##### Step 4: Check you have all the required executables
`seqenv` will search for several different binaries as it processes your data. Please check all of these are available in your `$PATH`:

    $ which blastn

##### Step 5: Get a local copy of the NT database
You can choose the database you want to BLAST against. By default we will search against `nt`. So check your `~/.ncbirc` file for the adequate references. If you don't have a copy of such a database, you need to download one to your machine.

##### Step 6: Install all required python packages
`seqenv` uses several third party python libraries. You can get them by running these commands:

    $ pip install biopython
    $ pip install sh
    $ pip install pandas
    $ pip install tqdm
    $ pip install biom-format
    $ pip install requests
    $ pip install pygraphviz
    $ pip install networkx
    $ pip install Orange-Bioinformatics

If you are on a machine that does not authorize you to install packages like that you can try to install them by adding the `--user` option.

If you are using a python manager such as pyenv, don't forget to rehash the binary links at the end:

    $ pyenv rehash

Now, you can check that it all works like this:

    $ python -c "import seqenv"

### Further documentation
The code is written in a clean, object-oriented, and pythonic way. It is also carries a fair deal of documentation. You can see this within the code itself or by going here: http://xapple.github.io/seqenv