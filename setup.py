#!/usr/bin/env python

'''
setup.py file for Scavager
'''
from setuptools import setup
#version = open('VERSION').readline().strip()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name             = 'Biosaur',
    version          = '0.0.1',
    description      = '''Proteomics post-search algorithm''',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    author           = 'Mark Ivanov & Lev Levitsky & Julia Bubis',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'https://bitbucket.org/markmipt/scavager',
    packages         = ['biosaur_src', ],
    install_requires = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 3',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    license          = 'License :: OSI Approved :: Apache Software License',
    entry_points     = {'console_scripts': ['biosaur = biosaur_src.biosaur:run']}
    )
