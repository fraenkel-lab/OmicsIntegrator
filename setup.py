#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    #requirements?
]

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
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    data_files=[('data',[]),
                 ('data/matrix_files',[]),
                 ('example',[]),
                 ('example/emt',[])],
    scripts=['bin/garnet.py','bin/motif_fsa_scores.py','bin/get_window_binding_matrix.py','bin/motif_regression.py','bin/map_peaks_to_known_genes.py','bin/zipTgms.py']
    )
