import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d,UnivariateSpline
from abc import ABC, abstractmethod

class ABCRate(ABC):
    range = (-np.inf, np.inf)
    @abstractmethod
    def integral(self, t0:float ,t1:float) -> float:
        return 0

    def total(self) -> float:
        return self.integral(*self.range)
        

class _limited(ABCRate):
    def __init__(self, r0, new_range):
        self._r0 = r0
        self.range = new_range
    def __call__(self, t):
        return (self.range[0]<=t)*(t<=self.range[1])*self._r0(t)
    def integral(self, t0,t1):
        R0,R1 = self.range
        t0 = min(max(t0,R0),R1)
        t1 = max(min(t1,R1),R0)
        return self._r0.integral(t0,t1)

class _mul(ABCRate):
    def __init__(self, r0, factor):
        self._r0 = r0
        self.C = factor
        self.range = r0.range
    def __call__(self, t):
        return self.C*self._r0(t)
    def integral(self, t0,t1):
        return self.C*self._r0.integral(t0,t1)

class _sum(ABCRate):
    def __init__(self, r0, r1):
        self.r0 = r0
        self.r1 = r1
        self.range = (min(r0.range[0], r1.range[0]),
                      max(r0.range[1], r1.range[1]))

    def __call__(self, t):
        return self.r0(t)+self.r1(t)
    def integral(self, t0,t1):
        return self.r0.integral(t0,t1)+self.r1.integral(t0,t1)
    
class _shift(ABCRate):
    def __init__(self, r0, dt):
        self.r0=r0
        self.dt=dt
        self.range = (r0.range[0]+dt,r0.range[1]+dt)
    def __call__(self, t):
        return self.r0(t-self.dt)
    def integral(self, t0,t1):
        return self.r0.integral(t0-self.dt,t1-self.dt)

class _invert(ABCRate):
    def __init__(self, r0):
        self.r0=r0
        self.range = (-r0.range[1],-r0.range[0])
    def __call__(self, t):
        return self.r0(-t)
    def integral(self, t0,t1):
        return self.r0.integral(-t1,-t0)

class Const(ABCRate):
    """Constant rate in time
    
    Args:
        c(float): constant rate value
    """
    def __init__(self, c):
        self.c=c
    def __call__(self, t):
        return self.c*np.ones_like(t)
    def integral(self, t0,t1):
        return self.c*(t1-t0)

class Func(ABCRate):
    """Rate defined by the function
    
    Args:
        f(callable[float]->float): the desired rate vs. time function
    """
    def __init__(self, f):
        self.f=f
    def __call__(self,t):
        return self.f(t)
    def integral(self, t0,t1):
        return quad(self.f,t0,t1)[0]

class Interpolated(ABCRate):
    """Rate defined by linear interpolation of the given points
    
    Args:
        x(1D array-like): time values
        y(1D array-like): rate values
        **kwargs(dict): additional parameters to pass to :class:`scipy.interp.UnivariateSpline`
            default = `dict(ext=1, k=1, s=1)` is for linear interpolation.
    """
    def __init__(self,x,y,**kwargs):
        self.range = (x[0],x[-1])
        kwargs.setdefault('ext',1)
        kwargs.setdefault('k',1)
        kwargs.setdefault('s',0)
        self.f = UnivariateSpline(x,y,**kwargs)
    def __call__(self,t):
        return self.f(t)
    def integral(self, t0,t1):
        return self.f.integral(t0,t1)

    
def _sort(x,y):
        idx=np.argsort(x)
        return x[idx],y[idx]
    
class LogRate(ABCRate):
    def __init__(self, x,y, extrapolate=False):
        x,y = _sort(x,y)
        def __calc_log_params(x,y):
            a = np.diff(np.log(y))/np.diff(np.log(x))
            yi = (x[1:]*y[1:]-y[:-1]*x[:-1])/(a+1)
            yi = np.nan_to_num(yi)
            yi=np.append([0],yi)
            yi = np.cumsum(yi)
            return a, yi

        self.x,self.y = x,y
        self.a,self.yi = __calc_log_params(x,y)
        self.extrapolate = extrapolate
        self.idx_max = len(self.a)-1
        self.range=(min(x),max(x))

    def _index(self,x):
        if np.isscalar(x):
            x = [x]
        idx = np.searchsorted(self.x,x)-1
        if self.extrapolate:
            idx[idx<0]=0
            idx[idx>self.idx_max] = self.idx_max
        return idx
    def __call__(self, x):
        idx = self._index(x)
        return self._eval(idx)(x)
    
    def _eval(self, idx):
        a = self.a[idx]
        x0,y0 = self.x[idx],self.y[idx]
        return lambda x:y0*(x/x0)**a
    
    def _int(self, x):
        idx = self._index(x)
        a = self.a[idx]
        x0,y0 = self.x[idx],self.y[idx]
        y1 = self._eval(idx)(x)
        res = self.yi[idx]+(y1*x-y0*x0)/(a+1)
        return np.squeeze(res)
    
    def integral(self,x0,x1):
        return self._int(x1)-self._int(x0)
ABCRate.__add__ = lambda self, other: _sum(self,other)
ABCRate.__mul__ = lambda self, factor:_mul(self,factor)
ABCRate.__rmul__= lambda self, factor:_mul(self,factor)
ABCRate.shift   = lambda self, dt: _shift(self,dt)
ABCRate.invert  = lambda self: _invert(self)

def rate(a, *, range=None):
    """ create a Rate object

    Args:
        a (`scalar` or `callable` or tuple[floats,floats] or :class:`ABCRate`):
            Object to use to create the rate
        range (`None` or `tuple(float,float)`):
            if given tuple (x0,x1), the rate outside (x0,x1) will be 0
    Returns:
        a constructed rate object (inherited from :class:`ABCRate`).
        Depending on the input type will produce following classes:

        * `scalar`   - constant rate :class:`sn_stat.rate.Const`
        * `callable` - rate defined by the function :class:`sn_stat.rate.Func`
        * `tuple(x,y)` - :class:`sn_stat.rate.Interpolated` (x,y)
        * :class:`ABCRate`: use given rate

    """
    def __make_rate(a):
        if isinstance(a,ABCRate):
            return a
        if callable(a):
            return Func(a)
        elif np.isscalar(a):
            return Const(a)
        else:
            return Interpolated(*a)

    r = __make_rate(a)
    if(range is not None):
        r = _limited(r, new_range=range)
    return r
    
def log_rate(a):
    """ Create a Rate object with log-log interpolation, 

    This is the alternative to :func:`rate` method for interpolation, 
    using the log-log interpolation instead of a linear one.
    I.e. assuming that between points (x0,y1) and (x1,y1)

    .. math:: y(x) = y_0(x/x_0)^{\\alpha}

    where :math:`\\alpha = \\ln(y_1/y_0)/\\ln(x_1/x_0)`
    
    Args:
        a: tuple of 1-D array-like (x,y)
            x values should be sorted and have the same sign (all x<0 or x>0)
    Returns:
        rate object
    Raises:
        ValueError: if not all x are of the same sign

    """
    x,y = a

    if np.all(x>=0):
        return LogRate(x,y)
    elif np.all(x<=0):
        return LogRate(-x,y).invert()
    else:
        raise ValueError("x values should have same sign")
    
