.. SN Stat documentation master file, created by
   sphinx-quickstart on Fri Feb 26 21:02:18 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SuperNova statistical analysis package (`sn_stat`)
==================================================

``sn_stat`` provides the implementation of the statistical methods used for search for the supernova neutrino signal in the series of neutrino interactions.


Task
----

Experiment measures the neutrino interactions in the detector, recording their timestamps :math:`\{t_i\}`.

The statistical analysis needs to determine if the measured interactions agree with the background rate, or if they contain a signal from the supernova, and, potentially, to identify the time of the signal start. 

When there is no supernova, the expected event rate at time *t* is :math:`B(t)` events/s.
In case of supernova starting at time :math:`t_0`, the event rate vs. time is :math:`B(t)+S(t-t_0)`.

The signal **significance** value is calculated for every assumed supernova time :math:`t_0`, measuring the deviation from the background-only hypothesis.

What ``sn_stat`` can do:
------------------------

* Define and manipulate the event rates vs. time :func:`sn_stat.rate`
* Calculate the Log Likelihood Ratio :class:`sn_stat.LLR`
* Calculate the significance of the supernova signal observation in the measured data :class:`sn_stat.ShapeAnalysis`,:class:`sn_stat.CountAnalysis`
* * And the expected distribution of the significance for a given hypothesis (rate)
* Generate the events timestamps based on the event rate :class:`sn_stat.Sampler`


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   usage
   reference


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
