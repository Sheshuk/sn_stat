from sn_stat import rate, LLR
from sn_stat.llr import Distr
from hypothesis import strategies as st, given, assume
from hypothesis.extra.numpy import arrays
import numpy as np
import pytest

Tvalue = st.floats(-10000,10000)
Trange = st.tuples(Tvalue, Tvalue).map(sorted)

@given(Tvalue,Tvalue, Trange)
def test_llr_const_single(t0,t_data,time_window):
    with pytest.raises(ValueError):
        llr = LLR(1,2)
    
    l = LLR(1,1,time_window=time_window)
    if(time_window[0] <= t_data-t0 <= time_window[1]):
        assert l(t_data,t0)==np.log(2)
    else:
        assert l(t_data,t0)==0

Xval = st.floats(-10000,10000)
Yval = st.floats(0,10000)

@st.composite
def distrS(draw,size=st.integers(1,100),bins=Xval,vals=Yval):
    n = draw(size)
    xs = draw(arrays(float,elements=bins, shape=n+1, unique=True))
    xs.sort()
    ys = draw(arrays(float,elements=vals, shape=n))
    assume(np.any(ys>0))
    ys = ys/ys.sum()
    d = Distr(xs,ys)
    d.set_interpolation()
    return d

@given(distrS())
def test_distr_sf(d):
    assert np.allclose(d.sf(-np.inf),1)
    assert np.allclose(d.sf( np.inf),0)
    
    assert np.allclose(d.sf(d.bins.min()),1)
    assert np.allclose(d.sf(d.bins.max()),0)
    
    assert np.allclose(d.sf(d.bins[1:]), 1.-np.cumsum(d.vals))

@given(distrS())
def test_distr_isf(d):
    assert d.isf(1) == d.bins.min()
    assert d.isf(0) == d.bins.max()

@given(distrS())
def test_distr_hist(d):
    vals,bins = d.histogram(d.bins)
    assert np.allclose(bins,d.bins)
    assert np.allclose(vals,d.vals)
