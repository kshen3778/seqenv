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
    long_description_content_type = 'text/markdown',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_packages(),

    include_package_data = True,
    scripts          = ['seqenv/seqenv'],
    install_requires = ['numpy==1.24.3', 'matplotlib==3.7.1', 'biopython==1.81', 'sh==2.0.4', 'pandas==2.0.2', 'tqdm==4.65.0', 'biom-format==2.1.15',
                        'requests==2.31.0', 'pygraphviz==1.11', 'networkx==3.1', 'Orange-Bioinformatics==2.6.25'],
)
