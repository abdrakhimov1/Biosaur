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
    version          = '2.0.3',
    description      = '''Proteomics post-search algorithm''',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    author           = 'Abdrakhimov Daniil & Mark Ivanov',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'https://github.com/abdrakhimov1/Biosaur',
    packages         = ['biosaur_src', ],
    install_requires = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 3',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    license          = 'License :: OSI Approved :: GNU General Public License',
    entry_points     = {'console_scripts': ['biosaur = biosaur_src.biosaur:run']}
    )
