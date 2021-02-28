#!/usr/bin/env python

from distutils.core import setup
import sn_stat

setup(name='sn-stat',
        version=sn_stat.__version__,
        description='Statistical methods for supernova neutrino detection',
        author='Andrey Sheshukov',
        author_email='ash@jinr.ru',
        licence='GNU GPLv3',
        packages=['sn_stat'],
        install_requires=['numpy','scipy'],
        extras_require={'doc':['sphinx','sphinx-rtd-theme'],
                        'test':['pytest','hypothesis']}
     )

