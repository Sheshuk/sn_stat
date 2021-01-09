import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

class Sampler:
    def __init__(self,r, time_window=[0,10], Npoints=1000):
        x = np.linspace(*time_window,Npoints)
        y = r(x)
        self._set(x,y)
    
    def _set(self,x,y):
        order = np.argsort(x)
        self.x=x[order]
        self.y=y[order]
        ycum = cumtrapz(y=self.y,x=self.x, initial=0)
        self.ycum=ycum
        Ytotal = ycum[-1]
        #print(Ytotal)
        p= ycum/Ytotal
        self.x_of_p = interp1d(p,self.x)
        self.Ytotal = Ytotal
        
    def sample(self):
        Ntot = np.random.poisson(self.Ytotal)
        ps = np.random.rand(Ntot)
        return self.x_of_p(ps)
    
