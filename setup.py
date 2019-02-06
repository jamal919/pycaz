# -*- coding: utf-8 -*-
"""
Setup for pyschism

@author: khan
"""

import setuptools

setuptools.setup(
    name='pyschism',
    version='0.1',
    author='Jamal Khan',
    author_email='jamal919@gmail.com',
    description='Python Package for Handling SCHISM Model',
    packages=setuptools.find_packages(),
    license='GPL v3',
    install_requires=[
        'numpy',
        'matplotlib',
        'basemap',
        'glob',
    ],
    include_package_data=True,
    url="https://github.com/jamal919/pyschism",
    long_description=open('README.md').read(),
)