#!/usr/bin/env python

from distutils.core import setup

setup(name='sn-stat',
        version='0.3.2',
        description='Statistical methods for supernova neutrino detection',
        author='Andrey Sheshukov',
        author_email='ash@jinr.ru',
        licence='GNU GPLv3',
        packages=['sn_stat'],
        install_requires=['numpy','scipy'],
        extras_require={'doc':['sphinx','sphinx-rtd-theme'],
                        'test':['pytest','hypothesis','flake8']}
     )

