# -*- coding: utf-8 -*-
"""
Setup for pycaz

@author: khan
"""

import setuptools

setuptools.setup(
    name='pycaz',
    version='0',
    author='Jamal Khan',
    author_email='4151009+jamal919@users.noreply.github.com',
    description='A python based analysis toolbox',
    packages=setuptools.find_packages(),
    license='GPL v3',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas'
    ],
    include_package_data=True,
    url="https://github.com/jamal919/pycaz",
    long_description=open('README.md').read(),
)