#!/usr/bin/env python

import setuptools 

with open('README.md') as f:
    readme = f.read()

setuptools.setup(name='sn-stat',
        description='Statistical methods for supernova neutrino detection',
        long_description=readme,
        long_description_content_type='text/markdown',
        author='Andrey Sheshukov',
        author_email='ash@jinr.ru',
        licence='GNU GPLv3',
        packages=['sn_stat'],
        install_requires=['numpy','scipy'],
        extras_require={'doc':['sphinx','sphinx-rtd-theme'],
                        'test':['pytest','hypothesis','flake8']},
        python_requires='>=3.7'
     )

