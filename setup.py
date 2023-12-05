# -*- coding: utf-8 -*-
"""
Setup for pycaz

@author: khan
"""

import setuptools

setuptools.setup(
    name='pycaz',
    author='Jamal Khan',
    author_email='4151009+jamal919@users.noreply.github.com',
    description='A python based analysis toolbox',
    packages=setuptools.find_packages(),
    license='GPL v3',
    setup_requires=['setuptools-git-versioning', 'setuptools_scm'],
    setuptools_git_versioning={
        "enabled": True,
        "template": "{tag}",
        "dev_template": "{tag}.post{ccount}+git.{sha}",
        "dirty_template": "{tag}.post{ccount}+git.{sha}.dirty",
        "starting_version": "0.1"
    },
    use_scm_version=True,
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'xarray',
        'requests',
        'beautifulsoup4',
        'utide'
    ],
    include_package_data=True,
    url="https://github.com/jamal919/pycaz",
    long_description=open('README.md').read()
)