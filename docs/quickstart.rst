Quickstart 
==========

Define the event rate
---------------------
.. code-block:: python

    from sn_stat import rate
    
    B = rate(0.5) #constant background rate
    S = rate(([0,1,10],[0,5,0]) #interpolated rate

For more information on defining rates see :ref:`event_rates`

Generate data sample
--------------------
:class:`sn_stat.Sampler` can produce random event samples from a given event rate

.. code-block:: python

    from sn_stat import Sampler
    tSN_true = 5
    smplr = Sampler(B+S.shift(5), time_window=[-10,20] ) #will generate events from signal starting at `t=5`     

Significance on a given data set
--------------------------------
.. code-block:: python

    from sn_stat 
