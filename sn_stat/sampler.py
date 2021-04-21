import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

class Sampler:
    """ Generates random event samples (timestamps) following the given event rate"""
    def __init__(self,r, time_window=[0,10], Npoints=1000):
        """
        Args:
            r (rate): event rate
            time_window (tuple[float,float]): limits in which the events are generated
            Npoints (int): number of subdivisions for the integration.
                The initial distribution is approximated by `Npoints` 
                equidistant points in the `time_window` range
        """
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
        """ Produce the random events

        Returns:
            ndarray: 1-d array with the events timestamps
        """
        Ntot = np.random.poisson(self.Ytotal)
        ps = np.random.rand(Ntot)
        return self.x_of_p(ps)
    
