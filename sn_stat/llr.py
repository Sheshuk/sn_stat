import numpy as np
from scipy import stats, fft, interpolate
from .rate import rate

class Distr:
    def __init__(self,bins,vals):
        self.bins=bins
        self.vals=vals
    def set_interpolation(self):
        def make_func(x,y):
            X,Y = x,y
            return interpolate.interp1d(X,Y, fill_value=(Y[0],Y[-1]), bounds_error=False, assume_sorted=True, kind='next')

        tail = np.cumsum(self.vals[::-1])[::-1]
        self.tail = np.concatenate((tail,[0]))
        self.vals = np.concatenate(([0],self.vals,[0]))
        self.binc = np.concatenate(([self.bins[0]],
                                     0.5*(self.bins[1:]+self.bins[:-1]),
                                     [self.bins[-1]]))
            
        #set functions
        self.pdf  = make_func(x=self.binc,y=self.vals)
        self.sf   = make_func(x=self.bins,y=self.tail)
        self.isf  = make_func(y=self.bins[::-1],x=self.tail[::-1])
    def histogram(self, bins):
        if np.isscalar(bins):
            bins = np.linspace(self.bins[0],self.bins[-1],bins)
        N = -np.diff(self.sf(bins))
        return N, bins
    
class LLR:
    """Log Likelihood Ratio(LLR) calculator

        Log likelihood ratio for H0 (B) and H1 (B+S) hypotheses:

            L(t,t0) = log(1+S(t-t0)/B(t))

        where t is the event time and t0 is assumed signal start time.

    """
    def __init__(self, S, B, time_window=None):
        """
        parameters:
        -----------
        S: rate|float|function|(x,y) tuple
            expected signal event rate vs. time
        B: rate|float|function|(x,y) tuple
            background rate vs. time
        time_window: tuple(T0,T1) or None
            the limits around t0 in which to take the signal.
            if None then try to take the full range from S (via S.range)

        Usually time_window should be the same as the range of S, or smaller if you want to consider only part of the signal shape.

        """
        self.S = rate(S)
        if(time_window is None):
            time_window = S.range

        if np.any(np.isinf(time_window)):
            raise ArgumentError(f'Cannot work with infinite time window: {time_window}')

        self.S0=self.S.integral(*time_window)
        self.B = rate(B) 

        self.time_window=np.array(time_window)
    
    def llr(self,ts,t0):
        if ts.size==0: 
            return np.zeros((1,len(t0)))
        tSN = ts-np.expand_dims(t0,1)
        res = np.log(1+self.S(tSN)/self.B(ts))
        res[(tSN<self.time_window[0])|(tSN>self.time_window[1])]=0
        return res
        
    def __call__(self,ts,t0):
        res = self.llr(ts,t0)
        return np.sum(res, axis=1)

    def sample(self,hypothesis, Npoints,t0):
        #sample the LLR with hypothesis
        ts = np.linspace(*self.time_window,Npoints)+t0
        ls = self.llr(ts,t0=[t0]).flatten()
        ws = hypothesis(ts)
        return ls,ws
    def l_range(self, t0, Npoints=10000):
        """
        return (min, max) LLR values for given t0
        """
        ls,_ = self.sample(hypothesis=self.B, Npoints=Npoints, t0=t0)
        return min(ls),max(ls)
    
    def distr(self, hypothesis='H0', t0=0, normal=False, Npoints=10000, l_bin_width=1e-3, **kwargs):
        if hypothesis=='H0':
            hypothesis=self.B
        ls,ws = self.sample(hypothesis,Npoints,t0)
        if normal:
            ws/=ws.sum()
            mu  = ls@ws
            var = (ls**2)@ws-mu**2
            var = max(var,1e-16)
            return stats.norm(loc=mu,scale=np.sqrt(var))
        # define the binning 
        binsl = np.arange(0,ls.max()+2*l_bin_width,l_bin_width)-l_bin_width/2.
        # produce the distribution
        H1, binl = np.histogram(ls, weights=ws, bins=binsl, density=False)
        H1/=H1.sum()
        return Distr(bins = binl, vals=H1)
    
def JointDistr(llrs, hypos='H0', t0=0, R_threshold=100, dl=1e-3, **kwargs):
    """
    Calculate the joint distribution of `llrs` under hypotheses `hypos`
    Parameters:
    * llrs -  an iterable of `LLR` objects
    * hypos - an iterable of `rate` objects
              or 'H0' string, taking background rates from each LLR
    * t0    - time of expected SN start (for LLR calc)
    * R_threshold - if the integrated rate in hypothesis is above this threshold,
                    a Gaussian approximation is used for this distribution.
                    If the rate is below - a precise calculation with FFT is performed.
    * dl    - LLR bin size for distributions. Ignored, if all distributions are gaussian.
    * **kwargs -arguments that will be passed to LLR.distr() method.
    
    Returns: Distr for the joint (sum) LLRs
    """
    def NormDistr(distrs,R,**kwargs):
        if len(distrs)==0:
            return None
        mu  = np.array([d.mean() for d in distrs])
        var = np.array ([d.var()  for d in distrs])
        return stats.norm(loc=mu@R,scale=np.sqrt((mu**2+var)@R))

    def FFTDistr(distrs, R, t0=0, dl=1e-3, epsilon=1e-16, **kwargs):
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

    def combine_distrs(*ds,Nbins=1000,epsilon=1e-16):
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
        hypos = [l.B for l in llrs]
    #prepare the rates for each experiment
    R = np.array([h.integral(*l.time_window+t0) for l,h in zip(llrs,hypos)])
   
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
 
    H1s=np.array([l.distr(h,t0,l_bin_width=llr_step,normal=is_norm,**kwargs) for l,h,is_norm in zip(llrs,hypos, largeR)])
    d1 = NormDistr(H1s[largeR], R[largeR],**kwargs)
    d2 = FFTDistr (H1s[smallR], R[smallR],dl=llr_step, **kwargs)
    return combine_distrs(d1,d2)
