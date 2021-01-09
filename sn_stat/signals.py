import numpy as np
from .rate import rate

class Signal:
    def __init__(self, s, distance=1):
        self.s = rate(s)
        self.d0 = distance
    def at(self,distance):
        return self.s*(self.d0/distance)**2

def ccSN(S0=1,t_rise=0.1,t_decay=1):
    norm = (t_decay**2)/(t_rise+t_decay)
    def _s(t):
        #res = np.zeros_like(t)
        #idx = np.where(t>0)
        #res[idx] = S0/norm*(1-np.exp(-t[idx]/t_rise))*np.exp(-t[idx]/t_decay)
        T = np.where(t>0,t,0)
        res = np.where(t>0, S0/norm*(1-np.exp(-T/t_rise))*np.exp(-T/t_decay), 0)
        #res = np.array(res)
        #res[res<0]=0
        return res
    return Signal(_s, distance=1)

def preSN(S0=1,t_rise=100,t_decay=0.1):
    f = ccSN(S0,t_rise=t_decay,t_decay=t_rise)
    return Signal(lambda t:f(t_decay-t), distance=1)

def from_file(fname, dt=5e-3, distance=10, scale=1):
    s = np.loadtxt(fname)*scale/dt
    t = np.arange(len(s))*dt
    return Signal((t,s), distance=distance)