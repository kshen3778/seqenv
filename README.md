# seqenv version 1.0.0
Assign environment ontology (EnvO) terms to short DNA sequences.

### Installing
To install `seqenv` onto your machine, use the python package manager:

    $ pip install seqenv

### Dependencies
* You need to have a copy nucleotide data base of NCBI (called `nt`) installed locally as well as the `blastn` executable in your `$PATH`.
* To compile parts of the software you will need a copy of the BOOST libraries as well as the SWIG libraries.
* The project also depends on some other python modules such as `biopython`. These will be installed automatically when calling the `pip` command.

### Usage
Once that is done, you can start processing FASTA files from the command line. For using the default parameters you can just type:

    $ seqenv sequences.fasta

We will then assume that you have inputed 16S sequences. To modify the database or use different type of sequence type:

    $ seqenv sequences.fasta --seqtype prot --db nr

To modify the minimum identity in the similarity search, use the following:

    $ seqenv sequences.fasta --identity 97

If you have abundance data you would like to add to your analysis you can specify it like this in a TSV file:

    $ seqenv sequences.fasta --abundances counts.tsv

### All parameters
   * `--seq_type`: Sequence type `nucl` or `prot` for nucleotides or amino acids, respectively (Default: `nucl`).
   * `--search_algo`: Search algorithm. Either `blast` or `usearch` (Default: `blast`).
   * `--search_db`: The database to search against (Default: `nt`). You can specify the full path or provide a `.ncbirc` file.
   * `--backtracking`: For every term identified by the tagger, we will propagate frequency counts up the acyclic directed graph described by the ontology. Defaults to `False`.
   * `--normalization`: Should we divide the counts of every input sequence by the number of text entries that were associated to it. Defaults to `True`.
   * `--num_threads`: Number of cores to use (Defaults to the total number of cores). Use 1 for non-parallel processing.
   * `--out_dir`: The output directory in which to store the result and intermediary files. Defaults to the same directory as the input file.
   * `--min_identity`: Minimum identity in similarity search (Default: 0.97). Note: not available when using `blastp`.
   * `--e_value`: Minimum e-value in similarity search (Default: 0.0001).
   * `--max_targets`: Maximum number of reference matches in similarity search (Default: 10).
   * `--min_coverage`: Minimum query coverage in similarity search (Default: 0.97).
   * `--abundances`: Abundances file (Default: None).
   * `--N`: If abundances are given, pick only the top N sequences (Default: 1000).

### Introduction
The continuous drop in the associated costs combined with the increased efficiency of the latest high-throughput sequencing technologies has resulted in an unprecedented growth in sequencing projects. Ongoing endeavours such as the [Earth Microbiome Project](http://www.earthmicrobiome.org) and the [Ocean Sampling Day](http://www.microb3.eu/osd) are transcending national boundaries and are attempting to characterise the global microbial taxonomic and functional diversity for the benefit of mankind. The collection of sequencing information generated by such efforts is vital to shed light on the ecological features and the processes characterising different ecosystems, yet, the full knowledge discovery potential can only be unleashed if the associated meta data is also exploited to extract hidden patterns. For example, with the majority of genomes submitted to NCBI, there is an associated PubMed publication and in some cases there is a GenBank field called "isolation sources" that contains rich environmental information.
With the advances in community-generated standards and the adherence to recommended annotation guidelines such as those of [MIxS](http://gensc.org/gc_wiki/index.php/MIxS) of the Genomics Standards Consortium, it is now feasible to support intelligent queries and automated inference on such text resources.

The [Environmental Ontology](http://environmentontology.org/) will be a critical part of this approach as it gives the ontology for the concise, controlled description of environments. It thus provides structured and controlled vocabulary for the unified meta data annotation, and also serves as a source for naming environmental information. Thus, we have developed the `seqenv` pipeline capable of annotating sequences with environment descriptive terms occurring within their records and/or in relevant literature. Given a set of sequences, `seqenv` retrieves highly similar sequences from public repositories (NCBI GenBank). Subsequently, from each of these records, text fields carrying environmental context information (such as the reference title and the **isolation source**) are extracted. Additionally, the associated **PubMed** links are followed and the relevant abstracts are collected. Once the relevant pieces of text for each matching sequence have been gathered, they are then processed by a text mining module capable of identifying EnvO terms mentioned in them. The identified EnvO terms along with their frequencies of occurrence are then subjected to clustering analysis and multivariate statistics. As a result, tagclouds and heatmaps of environment descriptive terms characterizing different sequences/samples are generated. The `seqenv` pipeline can be applied to any set of nucleotide and protein sequences. Annotation of metagenomic samples, in particular 16S rRNA sequences is also supported.

### Pipeline overview
![seqenv](https://bitbucket.org/repo/6g996b/images/3493861180-SEQenv.jpg "seqenv")

### Tutorial
We will first run `seqenv` on a 16S rRNA dataset using ***isolation sources*** as a text source. Here, `All_GoodT_C03.csv` is a species abundance file (3% OTUs) processed through [`AmpliconNoise`](https://code.google.com/p/ampliconnoise/) software and `All_GoodT_C03.fa` contains the corresponding sequences for the OTUs.

~~~
$ ls
All_GoodT_C03.csv
All_GoodT_C03.fa

$ seqenv -o 2 -n 1 -f All_GoodT_C03.fa -s All_GoodT_C03.csv -m 99 -q 99 -r 100
~~~

Once the pipeline has finished processing, you will have the following contents in the current folder:

~~~
$ ls
All_GoodT_C03_N1_blast_F_ENVO_OTUs.csv
All_GoodT_C03_N1_blast_F_ENVO_OTUs_labels.csv
All_GoodT_C03_N1_blast_F_ENVO_samples.csv
All_GoodT_C03_N1_blast_F_ENVO_samples_labels.csv
seqenv.log
~~~

-> #TODO

### Acknowledgments
`seqenv` was conceived and developed in the following hackathons supported by European Union's Earth System Science and Environmental Management ES1103 COST Action ("[Microbial ecology & the earth system: collaborating for insight and success with the new generation of sequencing tools](http://www.cost.eu/domains_actions/essem/Actions/ES1103)"):

- **From Signals to Environmentally Tagged Sequences** (Ref: ECOST-MEETING-ES1103-050912-018418), September 27th-29th 2012, Hellenic Centre for Marine Research, Crete, Greece.
- **From Signals to Environmentally Tagged Sequences II** (Ref: ECOST-MEETING-ES1103-100613-031037), June 10th-13th 2013, Hellenic Centre for Marine Research, Crete, Greece.
- **From Signals to Environmentally Tagged Sequences III** (Ref: XXX), September XXX 2014, Hellenic Centre for Marine Research, Crete, Greece.

This work would not have been possible without the advice and support of many people who attended the hackathons:

- [Umer Zeeshan Ijaz](http://userweb.eng.gla.ac.uk/umer.ijaz) (Umer.Ijaz@glasgow.ac.uk) [1,2]
- [Evangelos Pafilis](http://epafilis.info/) (pafilis@hcmr.gr) [1,2]
- [Chris Quince](http://www.gla.ac.uk/schools/engineering/staff/christopherquince/) (cq8u@udcf.gla.ac.uk) [2]
- Christina Pavloudi (cpavloud@hcmr.gr)
- Anastasis Oulas (oulas@hcmr.gr)
- Julia Schnetzer (jschnetz@mpi-bremen.de)
- Aaron Weimann (aaron.weimann@uni-duesseldorf.de)
- Alica Chronakova (alicach@upb.cas.cz)
- Ali Zeeshan Ijaz (alizeeshanijaz@gmail.com)
- Simon Berger (simon.berger@h-its.org)
- Lucas Sinclair (lucas.sinclair@me.com)

[1] Main developers
[2] Contact for correspondence

### News
* **August 2013**: Chris Quince presented a talk on `seqenv` at [STAMPS2013](https://stamps.mbl.edu/index.php/Main_Page). You can download the PDF of the presentation: [C Quince et. al., SeqEnv: Annotating sequences with environments (STAMPS 2013)](https://stamps.mbl.edu/images/4/44/Quince_SeqEnvSTAMPS2013.pdf)
