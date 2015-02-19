#!/usr/bin/env python
# -*- coding: utf-8 -*-
##This script collects all the files needed to build a wheel of garnet-forest:
##populate these data strucutres, then type 'setup.py bdist_wheel'

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = ["networkx","scipy","numpy","matplotlib"]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='garnet-forest',
    version='0.1.0',
    description="Python + Forest",
    long_description=readme + '\n\n' + history,
    author="Sara Gosline",
    author_email='sgosline@mit.edu',
    url='https://github.com/sgosline/garnetforest',
    packages=[
        'garnet-forest',
    ],
    package_dir={'garnet-forest':
                 'garnet-forest'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='garnet-forest',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Biologists with basic programming skills',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    data_files=[('data',['data/ucsc_hg19_kgXref.txt','data/ucsc_mm9_kgXref.txt',
                 'data/ucsc_hg19_knownGenes.txt','data/ucsc_mm9_knownGenes.txt','data/iref_mitab_miscore_2013_08_12_interactome.txt','data/README.rst']),
                 ('data/matrix_files',['data/matrix_files/vertebrates_clustered_motifs.txt','data/matrix_files/vertebrates_clustered_motifs_mIDs.txt','data/matrix_files/vertebrates_clustered_motifs_tfids.txt','data/matrix_files/vertebrates_clustered_motifs_up_tfids.txt','data/matrix_files/motif_thresholds.pkl','data/matrix_files/vertebrates_clustered_motifs.tamo']),
                 ('data/matrix_files/gifs',['data/matrix_files/gifs/*']),
                 ('example',['example/README.rst']),
                 ('example/a549',['example/a549/test-tgfb-data.py','example/a549/tgfb_forest.cfg',
'example/a549/Tgfb_exp.txt','example/a549/Tgfb_phos.txt','example/a549/tgfb_garnet.cfg','example/a549/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed','example/a549/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.fasta.gz']),
                 ('example/dnaseClus',['example/dnaseClus/dnaseClus_garnet.cfg','example/dnaseClus/wgEncodeRegDnaseClusteredV2.bed','example/dnaseClus/wgEncodeRegDnaseClusteredV2.fasta.gz']),
                 ('example/mcf7',['example/mcf7/mcf7_garnet.cfg','example/mcf7/wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.bed','example/mcf7/wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.fasta.gz']),
                 ('example/murineFib',['example/murineFib/murineFib_garnet.cfg','example/murineFib/wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.fasta.gz','example/murineFib/wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.narrowPeak'])],
    scripts=['scripts/garnet.py','scripts/motif_fsa_scores.py','scripts/get_window_binding_matrix.py','scripts/motif_regression.py','scripts/map_peaks_to_known_genes.py','scripts/zipTgms.py','scripts/forest.py']
    )
