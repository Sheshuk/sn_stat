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
    S_true=S.shift(tSN_true)
    R = B+S_true #rate with signal starting at `tSN_true`
    ts = Sampler(R, time_window=[-10,20]).sample() #generate events within given time range
    #or generate signal and bg events separately
    ts_B = Sampler(B,time_window=[-10,20]).sample() #generate the BG events!
    ts_S = Sampler(S_true,time_window=[-10,20]).sample() #generate the SG events!

Significance on a given data set
--------------------------------
.. code-block:: python

    from sn_stat 
