import numpy as np
from scipy import stats, fft, interpolate
from .det_config import DetConfig

class Distr:
    def __init__(self,bins,vals):
        self.bins=bins
        self.vals=vals
    def set_interpolation(self):
        def make_func(x,y):
            X,Y = x,y
            return interpolate.interp1d(X,Y, 
                    fill_value=(Y[0],Y[-1]), 
                    bounds_error=False, assume_sorted=True, kind='next')

        tail = np.cumsum(self.vals[::-1])[::-1]
        tail = np.concatenate((tail,[0]))
        vals = np.concatenate(([0],self.vals,[0]))
        binc = np.concatenate(([self.bins[0]],
                                     0.5*(self.bins[1:]+self.bins[:-1]),
                                     [self.bins[-1]]))
        #set functions
        self.pdf  = make_func(x=binc,y=vals)
        self.sf   = make_func(x=self.bins,y=tail)
        self.isf  = make_func(y=self.bins[::-1],x=tail[::-1])
    def histogram(self, bins):
        if np.isscalar(bins):
            bins = np.linspace(self.bins[0],self.bins[-1],bins)
        N = -np.diff(self.sf(bins))
        return N, bins
    def __repr__(self):
        return f'{__class__}(bins={self.bins}, vals={self.vals})'

class LLR:
    """ Log likelihood ratio for H0 (B) and H1 (B+S) hypotheses:

        .. math:: \\ell(t,t_0) = \\log\\left(1+\\frac{S(t-t_0)}{B(t)}\\right)

        where `t` is the event time and `t_0` is assumed signal start time.
        """
    def __init__(self, det: DetConfig):
        self.det = det
   
    def llr(self,ts,t0, w=1):
        if ts.size==0: 
            return np.zeros((1,len(t0)))
        tSN = ts-np.expand_dims(t0,1)
        res = np.log(1+self.det.S(tSN)/self.det.B(ts))*w
        res[(tSN<self.det.time_window[0])|(tSN>self.det.time_window[1])]=0
        return res
        
    def __call__(self,ts,t0, time_precision=None):
        """
        Calculate the LLR value for given set of measurements `ts`, assuming supernova times `t0`

        parameters
        ----------
        ts : iterable
            Measured interactions timestamps
        t0 : iterable
            Assumed supernova start times
        time_precision: float or `None`
            If not None: group the given `ts` to the time bins with given precision, 
            speeding up the calculation for large number of events

        returns
        -------
        llr : ndarray
            Cumulative LLR values for each value of `t0`

        """
        ts = np.array(ts, ndmin=1)
        t0 = np.array(t0, ndmin=1)
        if(time_precision):
            t,w = w,t = np.histogram(ts,bins=np.arange(ts.min(),ts.max(),time_precision))
            tc = 0.5*(t[1:]+t[:-1])
            tc = tc[w>0]
            w = w[w>0]
            res = self.llr(tc,t0,w)
        else:
            res = self.llr(ts,t0)
        return np.sum(res, axis=1)

    def sample(self,hypothesis, Nsamples,t0):
        #sample the LLR with hypothesis
        ts = np.linspace(*self.det.time_window,Nsamples)+t0
        ls = self.llr(ts,t0=[t0]).flatten()
        ws = hypothesis(ts)
        return ls,ws
    def l_range(self, t0, Nsamples=10000):
        """
        returns: (min, max) LLR values for given t0
        """
        ls,_ = self.sample(hypothesis=self.det.B, Nsamples=Nsamples, t0=t0)
        return min(ls),max(ls)
    
    def distr(self, hypothesis='H0', t0=0, *, normal=False, Nsamples=10000, dl='auto'):
        """
        Calculate the LLR distribution under given assumption of the event rate

        Args:
            hypothesis(`rate` or "H0"):
                the assumed event rate vs. time. If hypothesis=='H0' - use background rate (`self.det.B`)
            t0   (float): assumed supernova start time
            normal(bool): flag to use normal distribution, otherwise use precise FFT calculation

        Keyword Args:
            normal(bool):
                if True, approximate with normal distribution,
                otherwise calculate numerically with FFT
            dl(float or "auto"):
                LLR bin size for distribution.
                if 'auto' - set it to 1e-3*max(l)
                (ignored if `normal==True`)
            epsilon(float):
                calculation precision for FFT distributions
                (ignored if `normal==True`)
            Nsamples(int):
                number of points to sample the LLR values range
                (ignored if `normal==True`)

        Returns:
            distribution of the LLR values. 
            If `normal==True` returns :class:`scipy.stats.norm`,
            otherwise constructs :class:`Distr` with bins from 0 to maximum LLR value in given time window
        """
        if hypothesis=='H0':
            hypothesis=self.det.B
        ls,ws = self.sample(hypothesis,Nsamples,t0)
        if normal:
            ws/=ws.sum()
            mu  = ls@ws
            var = (ls**2)@ws-mu**2
            var = max(var,1e-16)
            return stats.norm(loc=mu,scale=np.sqrt(var))
        # define the binning 
        if dl=='auto':
            dl=ls.max()*1e-3
        binsl = np.arange(0,ls.max()+2*dl,dl)-dl/2.
        # produce the distribution
        H1, binl = np.histogram(ls, weights=ws, bins=binsl, density=False)
        H1/=H1.sum()
        return Distr(bins = binl, vals=H1)
    
def JointDistr(llrs, hypos='H0', t0=0, R_threshold=100, *, dl=1e-3, epsilon=1e-16, Nsamples=10000):
    """
    Calculate the joint distribution of `llrs` under hypotheses `hypos`

    Args:
        llrs (iterable of :class:`LLR`): 
            configurations for each experiment
        hypos (iterable of `rate` or  "H0"): 
            expected rate for each experiment,
            or "H0" string - taking background rates from each LLR
        t0 (float): 
            time of assumed SN start 
        R_threshold (int): 
            if the integrated rate in hypothesis is above this threshold,
            a Gaussian approximation is used for this distribution,
            otherwise a precise calculation with FFT is performed.

    Keyword Args:
        dl(float):
            LLR bin size for distributions. Ignored, if all distributions are gaussian.
        epsilon(float):
            calculation precision for FFT distributions
            (ignored if all distrs are gaussian)
        Nsamples(int):
            number of points to sample the LLR values range
            (ignored if all distrs are gaussian)
    
    Returns:
        :class:`Distr`: 
            a distribution for the joint (sum) of individual LLRs under the given hypothesis
    """
    def NormDistr(distrs,R):
        if len(distrs)==0:
            return None
        mu  = np.array([d.mean() for d in distrs])
        var = np.array ([d.var()  for d in distrs])
        return stats.norm(loc=mu@R,scale=np.sqrt((mu**2+var)@R))

    def FFTDistr(distrs, R, t0=0, dl=1e-3):
        if len(distrs)==0:
            return None
        #determine the resulting number of bins as a maximum of 
        Ns = 2*stats.poisson.isf(mu=R,q=epsilon)
        npoints = Ns*np.array([len(H1.vals)-1 for H1 in distrs])
        nbins = int(npoints.sum())+1

        # calculate fourier transform
        Hz = [fft.fft(H1.vals, n=nbins) for H1 in distrs]
        Hz = np.array(Hz)

        FF = np.exp(R@(np.array(Hz)-1))
        vals =  np.real(fft.ifft(FF))
        vals[vals<epsilon]=0
        
        res = Distr(vals = vals,
                    bins = (np.arange(nbins+1)-0.5)*dl)
        res.set_interpolation()
        return res

    def combine_distrs(*ds,Nbins=1000):
        #remove "None" histos
        ds = [d for d in ds if d is not None]
        if(len(ds)==1):
            return ds[0] #only one distr

        #get ranges
        l_min = np.array([d.isf([1-epsilon]) for d in ds])
        l_max = np.array([d.isf([epsilon]) for d in ds])
        #calc bin size
        dl = np.min(l_max-l_min)/Nbins

        #calculate histograms
        bs = []
        hs = []

        for l0,l1,d in zip(l_min,l_max,ds):
            bins = np.arange(l0,l1,dl)
            h = -np.diff(d.sf(bins))
            bs+=[bins]
            hs+=[h]

        #convolution
        h = np.convolve(*hs)
        b = np.arange(len(h)+1)*dl + np.sum(l_min)
        res = Distr(b,h)
        res.set_interpolation()
        return res
    
    if hypos=='H0':
        hypos = [l.det.B for l in llrs]
    #prepare the rates for each experiment
    R = np.array([h.integral(*l.det.time_window+t0) for l,h in zip(llrs,hypos)])
   
    #divide small and large R cases
    largeR = (R>=R_threshold)
    smallR = largeR==False
    #prepare the distributions for each experiment)
    if np.any(smallR):
        llr_max = [l.l_range(t0=t0)[1] for l in np.array(llrs)[smallR]]
        llr_max = max(llr_max)
        llr_step=dl*llr_max
    else:
        llr_step=dl
 
    H1s=np.array([l.distr(h,t0,dl=llr_step,normal=is_norm,Nsamples=Nsamples) for l,h,is_norm in zip(llrs,hypos, largeR)])
    d1 = NormDistr(H1s[largeR], R[largeR])
    d2 = FFTDistr (H1s[smallR], R[smallR],dl=llr_step)
    return combine_distrs(d1,d2)

