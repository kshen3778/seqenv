from setuptools import setup, find_packages

setup(
    name             = 'seqenv',
    version          = '1.3.0',
    description      = 'Assign environment ontology (EnvO) terms to DNA sequences',
    license          = 'MIT',
    url              = 'https://github.com/xapple/seqenv',
    download_url     = 'https://github.com/xapple/seqenv/tarball/1.2.9',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    long_description = open('README.md').read(),
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_packages(),

    include_package_data = True,
    scripts          = ['seqenv/seqenv'],
    install_requires = ['numpy', 'matplotlib', 'biopython', 'sh', 'pandas', 'tqdm', 'biom-format',
                        'requests', 'pygraphviz', 'networkx', 'Orange-Bioinformatics'],
)
