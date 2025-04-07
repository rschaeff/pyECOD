#!/usr/bin/env python3
"""
Setup script for pyECOD
"""

from setuptools import setup, find_packages

setup(
    name="pyecod",
    version="0.1.0",
    description="Python toolkit for Evolutionary Classification of Protein Domains",
    author="RD Schaeffer",
    author_email="dustin.schaeffer@gmail.com",
    packages=find_packages(),
    install_requires=[
        "psycopg2-binary>=2.9.3",
        "pyyaml>=6.0",
        "biopython>=1.79",
        "numpy>=1.22.0",
        "pandas>=1.4.0",
        "matplotlib>=3.5.0",
        "requests>=2.27.0",
    ],
    scripts=[
        "ecod/scripts/manage_jobs.py",
    ],
    entry_points={
        'console_scripts': [
            'ecod=ecod.cli.main:main',
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
)