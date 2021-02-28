from scipy import stats
import numpy as np
from .llr import JointDistr

def p2z(p):
    "convert p-value to significance"
    return stats.norm.isf(p)
def z2p(z):
    "convert significance to p-value"
    return stats.norm.sf(z)

class ShapeAnalysis:
    """ Calculating significance using the shape analysis method.

    """
    def __init__(self, llrs, **params):
        self.llrs = llrs
        self.params=params
        self.d0 = JointDistr(self.llrs,**self.params)
    def __call__(self, data, t0):
        """
        calculate significance for the set of measurements
        """
        assert len(data)==len(self.llrs)
        ls = np.stack([l(d,t0) for l,d in zip(self.llrs, data)])
        return self.l2z(np.sum(ls, axis=0))
    
    def l_distr(self,hypos,add_bg=False):
        """
        calculate LLR distribution for given hypotheses"
        """
        assert len(hypos)==len(self.llrs)
        if(add_bg):
            hypos = [h+l.B for h,l in zip(hypos,self.llrs)]
        return JointDistr(self.llrs,hypos,**self.params)
    
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
    
    def l2p(self, l):
        "convert LLR to p-value"
        return self.d0.sf(l)
    def p2l(self, p):
        "convert p-value to LLR"
        return self.d0.isf(p)
    def l2z(self, l):
        "convert LLR to significance"
        return p2z(self.l2p(l))
    def z2l(self, z):
        "convert significance to LLR"
        return self.d0.isf(z2p(z))
    
