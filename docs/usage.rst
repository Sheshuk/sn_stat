Usage
=====

.. toctree::

    rate

Rate definition
---------------
For the shape analysis each experiment must define the background rate `B` 
and expected signal rate `S`. This is done with :func:`sn_stat.rate` function:

.. code-block:: python

   from sn_stat import rate
   B = rate(1) #constant rate 1 event/s
   
   #constant rate within given time range, 0 outside
   S1 = rate(2,range=[0,10]) 

   times = [0, 1., 10.]
   rates = [0, 10.,0.]
   S2 = rate((times,rates)) #triangular signal shape
   S3 = rate((times,rates), range=[0,5]) #triangular signal shape with reduced range



Detector configuration
----------------------

The :class:`st_stat.DetConfig` defines the configuration of the detector: background conditions, analysis time window, expected signal shape.

.. code-block:: python

   import sn_stat as st
   det1 = st.DetConfig(B=B, S=S1, time_window=[0,10])
   det2 = st.DetConfig(B=st.rate(10.), S=S2, time_window=[0,5])

Significance calculation
------------------------

There are currently two methods of significance calculation implemented: :class:`sn_stat.CountAnalysis` and :class:`sn_stat.ShapeAnalysis`.
They can be instantiated for the given detector configuration:

By default it can be used as 

.. code-block:: python

   ana1 = st.CountAnalysis(det1)
   ana2 = st.ShapeAnalysis(det2)

They can perform two tasks:

Process measured data 
*********************
To obtain test statistic value :meth:`sn_stat.Analysis.l_val` or supernova significance :meth:`sn_stat.Analysis.__call__`:

.. code-block:: python

   ts = ... #get events from the data
   t0s = np.arange(-10,10,0.01) #the assumed supernova positions we're interested in
   zs = ana1(ts,t0s)
   ls = ana1.l_val(ts,t0s)

Calculate expected distributions
********************************
We can calculate the expected distribution of the test statistics under given hypothesis
using :meth:`sn_stat.ShapeAnalysis.l_distr`, distribution of significance values using :meth:`sn_stat.ShapeAnalysis.z_distr`, and values, corresponding to certain probability quantiles :meth:`sn_stat.ShapeAnalysis.z_quant`.


