Usage
=====

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


Significance calculation
------------------------

The class :class:`sn_stat.ShapeAnalysis` is meant for the calculation of the supernova observation significance. 

By default it can be used as 

.. code-block:: python

   from sn_stat import LLR,ShapeAnalysis
   l1 = LLR(S=S1,B=1)
   l2 = LLR(S=S2,B=1)

   ana = ShapeAnalysis([l1,l2])


It can perform two kinds of calculation:

Significance quantiles
**********************
We can calculate the significance values, corresponding to certain significance quantiles, 
using :meth:`sn_stat.ShapeAnalysis.z_quant`:




Significance distribution
*************************
Method :meth:`sn_stat.ShapeAnalysis.z_distr`


Given set of measurements `ts`, calculate the 
    Calculate the significance distribution 
