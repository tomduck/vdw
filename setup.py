#! /usr/bin/env python3

import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

DESCRIPTION = """\
van der Waals fluid functions.
"""

VERSION = '0.1'

setup(name='vdw',
    version=VERSION,

    author='Thomas J. Duck',
    author_email='tomduck@tomduck.ca',
    description=DESCRIPTION,
    long_description=DESCRIPTION,
    license='GPL',
    keywords='van der Waals fluid functions',
    url='https://github.com/tomduck/vdw',

    install_requires=['numpy', 'matplotlib'],

    py_modules=['vdw'],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python'
        ]
    )
