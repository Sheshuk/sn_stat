.. _event_rates:

Event rates
===========

The event rate vs. time in most cases is defined with :func:`sn_stat.rate` function::

    from sn_stat import rate

Creation
--------

Rates can be constant in time::

    R0 = rate(1) #a constant rate of 1 evt/s

or as linear interpolation of the points::

    R1 = rate(([-1,0,10],[0,5,0])) #triangular signal shape

or from a function::

    R2 = rate(lambda x: np.sin(x/np.pi)**2)

Optionally you can provide `range` parameter, to limit the rate in time::

    R0_lim = rate(1, range=(0,1)) #will be 0 for x outside of the range

also you can take existing rates and limit them::

    R1_lim = rate(R1, range=(0,5)) #same as R1, but cropped

Operations
----------

Rates can be multiplied by factor, added together::

    B = rate(1) #background
    S = rate(10, range=(0,1)) #signal at certain distance
    R = B+S*4. #total event rate with the scaled signal

and shifted in time::

    S_shifted = S.shift(5) #Will start 5 seconds later


Then, rates can be used to calculate definite integrals::

    R0.integral(-10,10) #20 events
    R0_lim.integral(-10,10) #1 event

or get rate at the given time::
    
    R1(0)    #5 evt/s
    R1(-0.5) #2.5 evt/s
    R1(10)   #0 evt/s


Log-Log interpolation
---------------------

In some cases (like presupernova neutrino signal) it is feasible to use not the linear interpolation as in :class:`sn_stat.rate.Interpolated`, but a log-log interpolation. See :func:`sn_stat.log_rate`
