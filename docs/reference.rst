Reference
=========

rate
------------

.. autoclass:: sn_stat.rate
.. autoclass:: sn_stat.log_rate

Concrete rates
**************
.. automodule:: sn_stat.rate
    :members: Const,Func,Interpolated

sampler
-------
.. autoclass:: sn_stat.Sampler
    :members:

detector configuration
----------------------
.. autoclass:: sn_stat.DetConfig
    :members:

llr
-----------
.. automodule:: sn_stat.llr
    :members: JointDistr

.. autoclass:: sn_stat.llr.Distr
    :members:

.. autoclass:: sn_stat.LLR
    :special-members: __call__
    :members:

Analyses
--------------
.. autoclass:: sn_stat.CountingAnalysis
    :special-members: __call__
    :members:
    :inherited-members:


.. autoclass:: sn_stat.ShapeAnalysis
    :special-members: __call__
    :members:
    :inherited-members:
