#!/usr/bin/env python
# -*- coding: utf-8 -*-
##This script collects all the files needed to build a distribution of OmicsIntegrator:
##populate these data structures, then type 'setup.py sdist'

import sys, os

from setuptools import setup
from setuptools.command.test import test as TestCommand


with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = ["networkx","scipy","numpy","matplotlib"]

test_requirements = ['pytest']

gifdir='data/matrix_files/gifs/'
gif_files=[gifdir+g for g in os.listdir(gifdir)]

# Following pytest good practices
# https://pytest.org/latest/goodpractises.html
class PyTest(TestCommand):
    # Needed to set the --msgsteiner argument
    user_options = [('pytest-args=', 'a', '"Arguments to pass to py.test')]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ['tests']

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Import here, outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

setup(
    name='OmicsIntegrator',
    version='0.2.0',
    description="Python tool for high throughput data integration",
    long_description=readme + '\n\n' + history,
    author="Sara Gosline",
    author_email='sgosline@mit.edu',
    url='https://github.com/fraenkel-lab/OmicsIntegrator',
    packages=[],
#    package_dir={'./','./examples'},#'OmicsIntegrator':
                 #'OmicsIntegrator'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD License",
    zip_safe=True,
    keywords='OmicsIntegrator',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Biologists with basic programming skills',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    data_files=[('data',['data/ucsc_hg19_kgXref.txt','data/ucsc_mm9_kgXref.txt',
                 'data/ucsc_hg19_knownGenes.txt','data/ucsc_mm9_knownGenes.txt','data/iref_mitab_miscore_2013_08_12_interactome.txt','data/README.md']),
                 ('data/matrix_files',['data/matrix_files/vertebrates_clustered_motifs.txt','data/matrix_files/vertebrates_clustered_motifs_mIDs.txt','data/matrix_files/vertebrates_clustered_motifs_tfids.txt','data/matrix_files/vertebrates_clustered_motifs_up_tfids.txt','data/matrix_files/motif_thresholds.pkl','data/matrix_files/vertebrates_clustered_motifs.tamo']),
                 ('data/matrix_files/gifs',gif_files),
                 ('example',['example/README.md'])#,
#                 ('example/a549',['example/a549/test-tgfb-data.py','example/a549/tgfb_forest.cfg',
#'example/a549/Tgfb_exp.txt','example/a549/Tgfb_phos.txt','example/a549/tgfb_garnet.cfg','example/a549/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed','example/a549/wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.fasta.gz']),
#                 ('example/dnaseClus',['example/dnaseClus/dnaseClus_garnet.cfg','example/dnaseClus/wgEncodeRegDnaseClusteredV2.bed','example/dnaseClus/wgEncodeRegDnaseClusteredV2.fasta.gz']),
#                 ('example/mcf7',['example/mcf7/mcf7_garnet.cfg','example/mcf7/wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.bed','example/mcf7/wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.fasta.gz']),
#                 ('example/murineFib',['example/murineFib/murineFib_garnet.cfg','example/murineFib/wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.fasta.gz','example/murineFib/wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.narrowPeak'])
                 ],
    scripts=['scripts/garnet.py','scripts/motif_fsa_scores.py','scripts/get_window_binding_matrix.py','scripts/motif_regression.py','scripts/map_peaks_to_known_genes.py','scripts/zipTgms.py','scripts/forest.py'],
    tests_require=test_requirements,
    cmdclass={'test':PyTest}
)
