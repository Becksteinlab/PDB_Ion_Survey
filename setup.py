#! /usr/bin/python
"""Setuptools-based setup script for pdbionsurvey.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(name='pdbionsurvey',
      version='0.2.0-dev',
      description='survey of ion coordination in the Protein Data Bank',
      author='Kacey Clark',
      author_email='kacey.reidy@gmail.com',
      classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      scripts=[],
      license='GPLv3',
      long_description=open('README.md').read(),
      install_requires=['numpy', 'pandas', 'MDAnalysis', 'datreant', 'peakutils', 'matplotlib', 'seaborn', 'json', 'scipy', 'mmtf']
      )
