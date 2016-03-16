from distutils.core import setup

setup(
      name             = 'seqenv',
      version          = '1.0.7',
      description      = 'Assign environment ontology (EnvO) terms to short DNA sequences',
      license          = 'MIT',
      url              = 'https://github.com/xapple/seqenv',
      download_url     = 'https://github.com/xapple/seqenv/tarball/1.0.7',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['seqenv', 'seqenv.common', 'seqenv.fasta', 'seqenv.seqsearch'],
      scripts          = ['seqenv/seqenv'],
      install_requires = ['biopython', 'sh', 'pandas', 'tqdm', 'biom-format', 'requests'],
      long_description = open('README.md').read(),
      long_description_markdown_filename='README.md'
    )
