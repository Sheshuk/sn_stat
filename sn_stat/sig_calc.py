from scipy import stats
import numpy as np
from .llr import JointDistr, LLR
from . import DetConfig
from abc import ABC, abstractmethod
from scipy.stats import poisson
from collections.abc import Iterable

def p2z(p):
    "convert p-value to significance"
    return stats.norm.isf(p)
def z2p(z):
    "convert significance to p-value"
    return stats.norm.sf(z)


class Analysis(ABC):
    def __init__(self, discrete=False):
        self.d0 = self.l_distr(hypos="H0")
        if discrete:
            self._pmf = self.d0.pmf
        else:
            self._pmf = self.d0.pdf#lambda x:0

    @abstractmethod
    def l_distr(self, hypos, add_bg=False):
        """
        Calculate test statistic distribution

        Args:
            hypos: :class:rate or "H0"
                Hypothetical event we use to calculate the distribution
                If "H0" then use the background rate (self.det.B)

            add_bg: bool
                If set, then add the background rate to the `hypos`
        Returns:
            frozen poisson distribution
        """
        pass

    @abstractmethod
    def l_val(self, data, d0, **params):
        """
        Calculate test statistics value

        Args:
            data (iterable of array of float): 
                List with arrays of measured events time stamps for each detector.
                If there is only one detector, just an array(float) is enough
            t0 (ndarray of float):
                assumed time/times of signal start
            params:
                additional optional parameters for the current method
        Returns:
            ndarray of float:
                test statistic values for each value in `t0`
        """
        pass

    def __call__(self, data, t0, **params):
        """
        calculate significance for the set of measurements
        """
        return self.l2z(self.l_val(data,t0, **params))
 
    def l2p(self, l):
        "convert TestStatistics to p-value"
        return self.d0.sf(l)+self._pmf(l)
    def p2l(self, p):
        "convert p-value to TestStatistics"
        return self.d0.isf(p)
    def l2z(self, l):
        "convert TestStatistics to significance"
        return p2z(self.l2p(l))
    def z2l(self, z):
        "convert significance to TestStatistics"
        return self.d0.isf(z2p(z))

    def z_quant(self,hypos,add_bg=False,qs=z2p([-1,0,1])):
        """
        calculate zs, corresponding to quantiles of given hypotheses
        """
        d1 = self.l_distr(hypos,add_bg)
        ls = d1.isf(qs)
        zs = self.l2z(ls)
        return np.nan_to_num(zs)
    def z_distr(self,hypos, add_bg=False,zbins=100):
        """
        calculate significance distribution based on given hypotheses `hypos`
        """
        if (np.isscalar(zbins)):
            zbins=np.linspace(-5,5,zbins)
        lbins = self.z2l(zbins)
        d1 = self.l_distr(hypos,add_bg)
        
        return -np.diff(d1.sf(lbins)),zbins
    


class CountingAnalysis(Analysis):
    def __init__(self, det: DetConfig):
        """ 
        Calculating significance using the counting analysis method:
        using number of interactions within the time window (:class:`DetConfig.time_window`) as the test statistics (TS).

        Args:
            det(:class:`DetConfig`): 
                configuration for the experiment
        """
 
        self.det = det
        self.time_window = det.time_window
        super().__init__(discrete=True)

    def l_distr(self, hypos, add_bg=False):
        if(hypos=="H0"):
            hypos=self.det.B
        if(add_bg):
            hypos = hypos+self.det.B

        N = hypos.integral(*self.det.time_window)
        return poisson(mu=N)

    def l_val(self, data, t0, **params):
        data = np.array(data, ndmin=2).T
        t0 = np.array(t0, ndmin=2)
        tw = self.det.time_window
        T0,T1 = tw[0]+t0, tw[1]+t0
        return np.sum( (data>=T0)&(data<=T1), axis=0)


class ShapeAnalysis(Analysis):
    def __init__(self, detectors, **params):
        """ 
        Calculating significance using the shape analysis method, 
        using Log Likelihood Ratio (:class:`LLR`) 

        Args:
            detectors (single :class:`DetConfig` or iterable of :class:`DetConfig`): 
                configurations for each experiment. 
                Passing single DetConfig :code:`ShapeAnalysis(det)` is equivalent 
                to passing a list with one item :code:`ShapeAnalysis([det])`
                    
        Keyword Args:
            params (dict of kwargs):
                configuration arguments to be passed to :func:`sn_stat.llr.JointDistr`:
                    
        """
        if isinstance(detectors, DetConfig):
            detectors = [detectors]
        self.time_window = [
                min([d.time_window[0] for d in detectors]),
                max([d.time_window[1] for d in detectors])
                ]
        self.llrs = [LLR(d) for d in detectors]
        self.params=params
        self.det = detectors
        super().__init__()
    
    def l_val(self, data, t0, **params):
        if(len(data)!=len(self.llrs)):
            data = np.array(data, ndmin=2)
            assert data.shape[0]==len(self.llrs)
        ls = np.stack([l(d,t0,**params) for l,d in zip(self.llrs, data)])
        return np.sum(ls,axis=0)
   
    def l_distr(self,hypos,add_bg=False):
        if hypos!="H0":
            if not isinstance(hypos, Iterable):
                hypos=[hypos]
            assert len(hypos)==len(self.llrs)
            if(add_bg):
                hypos = [h+l.det.B for h,l in zip(hypos,self.llrs)]
        return JointDistr(self.llrs,hypos,**self.params)
    
