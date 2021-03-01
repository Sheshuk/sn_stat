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

Define LLR
----------
Background, expected signal and analysis time window should be defined for each considered experiment/detector/analysis.
These parameters should be stored in the Log Likelihood Ratio calculator class :class:`sn_stat.LLR`

.. code-block:: python

    from sn_stat import LLR

    #construct with existing rates
    B = rate(0.5) #constant background rate
    S = rate(([0,1,2],[0,5,0]) #interpolated rate
    llr = LLR(S, B) #define LLR calculator

    #alternatively rates can be constructed from given arguments on the fly
    llr = LLR(S=([0,1,2],[0,5,0]), B=0.5) #equivalent to previous

Significance on a given data set
--------------------------------
.. code-block:: python

    from sn_stat import ShapeAnalysis

    #analysis for a single experiment
    ana = ShapeAnalysis([llr]) #we pass one LLR calculator
    
    t0s = np.linspace(0,10) #assumed SN start time
    #Calculate the significance
    zs = ana(ts,t0s)

Significance quantiles and distribution
---------------------------------------
.. code-block:: python

    from sn_stat import rate, LLR, ShapeAnalysis
    
    B = rate(10)
    S = rate(([0,1,10],[0,10,0]))

    llr = LLR(S,B)
    ana = ShapeAnalysis([llr]) 
    
    #Calculate the significance quantiles
    from sn_stat import z2p
    quantiles =  z2p([-1,0,1]) #quantiles corresponding to median and +-sigma band

    zs0 = ana.z_quant([B], qs = quantiles) #calculate assumin only background
    print(zs0) #[-1,0,1] - zero significance in case of no supernova
    zs1 = ana.z_quant([B+S], qs = quantiles) #calculate assumin only background
    print(zs1) #[4.44822334, 5.76116125, 7.07409916] - high significance if the SN signal is seen


