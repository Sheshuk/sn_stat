#!/usr/bin/env python

from distutils.core import setup
from sn_stat import __version__

setup(name='sn-stat',
        version=__version__,
        description='Statistical methods for supernova neutrino detection',
        author='Andrey Sheshukov',
        author_email='ash@jinr.ru',
        licence='GNU GPLv3',
        packages=['sn_stat'],
        install_requires=['numpy','scipy'],
        extras_require={'doc':['sphinx','sphinx-rtd-theme'],
                        'test':['pytest','hypothesis']}
     )

