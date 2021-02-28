Quickstart 
==========

Define the event rate
---------------------
.. code-block:: python

    from sn_stat import rate
    
    B = rate(0.5) #constant background rate
    S = rate(([0,1,5],[0,5,0]) #interpolated rate

For more information on defining rates see :ref:`event_rates`

Generate data sample
--------------------
:class:`sn_stat.Sampler` can produce random event samples from a given event rate

.. code-block:: python

    from sn_stat import Sampler

    tSN_true = 5
    S_true=S.shift(tSN_true) #true signal starts at 5s
    ts = Sampler(B+S_true).sample() #generate events within given time range
    #or generate signal and bg events separately
    ts_B = Sampler(B).sample() #generate the BG events!
    ts_S = Sampler(S_true).sample() #generate the SG events!

Significance on a given data set
--------------------------------
.. code-block:: python

    from sn_stat 
