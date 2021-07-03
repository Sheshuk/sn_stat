from sn_stat import rate, log_rate
from hypothesis import strategies as st, given
from hypothesis.extra.numpy import arrays
import numpy as np

Xvalues = st.floats(-100000,100000)
Yvalues = st.floats(0,100000)

Ranges = st.tuples(Xvalues,Xvalues).map(sorted)

@given(c=Yvalues,xs=st.tuples(Xvalues,Xvalues))
def test_const(c,xs):
    r = rate(c)
    for x in xs:
        assert r(x)==c
    assert r.range == (-np.inf, np.inf)
    assert r.integral(xs[0],xs[1]) == c*(xs[1]-xs[0])

@given(c1=Yvalues, c2=Yvalues, xs=arrays(np.float64, elements=Xvalues, shape=(100)) )
def test_const_sum(c1,c2, xs):
    r1 = rate(c1)
    r2 = rate(c2)
    print(c1,c2)
    r = r1+r2
    assert np.allclose(r(xs),c1+c2)

# test ranges manipulation
@given(lims = Ranges, xs=st.lists(Xvalues))
def test_range(lims, xs):
    r = rate(1, range=lims)
    for x in xs:
        assert r(x) == (lims[0]<=x<=lims[1])


@given(Ranges, Ranges)
def test_range_add(rs0,rs1):
    r0 = rate(1, range=rs0)
    r1 = rate(2, range=rs1)
    assert r0.range==rs0
    assert r1.range==rs1
    r = r0+r1
    assert r.range[0] == min(rs0[0],rs1[0])
    assert r.range[1] == max(rs0[1],rs1[1])

@given(Ranges, Yvalues)
def test_range_mul(lims,factor):
    r0 = rate(1, range=lims)
    r1 = r0*factor
    assert r0.range==r1.range

@given(Ranges, Xvalues)
def test_range_shift(lims,dx):
    r0 = rate(1, range=lims)
    r1 = r0.shift(dx)
    assert r1.range[0]==r0.range[0]+dx
    assert r1.range[1]==r0.range[1]+dx

@given(lims=Ranges, t0t1=Ranges)
def test_range_integral(lims,t0t1):
    r = rate(1,range=lims)
    t0,t1 = t0t1
    if(lims[0]<=t0<=t1<=lims[1]):
        assert r.integral(t0,t1) == t1-t0
    elif(t0<=lims[0]<=t1<=lims[1]):
        assert r.integral(t0,t1) == t1-lims[0]
    elif(lims[0]<=t0<=lims[1]<=t1):
        assert r.integral(t0,t1) == lims[1]-t0
    elif(t0<=lims[0]<=lims[1]<=t1):
        assert r.integral(t0,t1) == lims[1]-lims[0]
    else:
        assert r.integral(t0,t1) == 0 

@st.composite
def xyS(draw, Xs=Xvalues, Ys=Yvalues, size=st.integers(2,100)):
    n = draw(size)
    x = draw(arrays(float,elements=Xs,shape=n, unique=True))
    x.sort()
    y = draw(arrays(float,elements=Ys,shape=n))
    return (x,y)

@given(xy=xyS())
def test_interp(xy):
    x,y = xy
    r = rate((x,y))
    assert r.range[0]==x[0]
    assert r.range[1]==x[-1]
    assert np.allclose(r(x), y)

