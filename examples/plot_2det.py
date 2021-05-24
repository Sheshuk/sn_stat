#!/bin/env python 

import sn_stat as sn
import numpy as np
import pylab as plt

#define the rates
sg1 = sn.signals.ccSN(S0=2)
sg2 = sn.signals.ccSN(S0=5, t_rise=2, t_decay=3)
#sg2 = sn.rate(([-5,0,1,5],[0,1,4,0]))

det1 = sn.DetConfig(B=0.1,S=sg1.at(1), time_window=[-5,10])
det2 = sn.DetConfig(B=1,S=sg2.at(1), time_window=[-5,10])

detectors = [det1,det2]
Tsn = 0.
T0,T1 = -60,60
Rs = []
for det in [det1,det2]:
    print(det)
    #total event rate vs time
    R = det.B+det.S.shift(Tsn)
    det.rate = R
    #create the dataset via sampler
    s = sn.Sampler(R, time_window=[T0,T1])
    det.data = s.sample()

#maximum time window:
tw = [  max(d.time_window[0] for d in detectors),
        min(d.time_window[1] for d in detectors)]

#assumed supernova start times
t0s = np.linspace(T0-tw[0], T1-tw[1],501)

for det in detectors:
    det.zs_s = sn.ShapeAnalysis   (det)(det.data,t0s)
    det.zs_c = sn.CountingAnalysis(det)(det.data,t0s)

data = [det.data for det in detectors]
zs_total = sn.ShapeAnalysis(detectors)(data,t0s)
#plot everything
plt.style.use('seaborn-talk')

nplots = 3
fig, ax = plt.subplots(nplots,1, figsize=(12,4*nplots), sharex=True, constrained_layout=True)

plt.sca(ax[0])
plt.ylabel('Event rate, 1/s')

for det in detectors:
    plt.sca(ax[0])
    plt.hist(det.data,bins=np.arange(T0,T1,1), density=False, histtype='step', lw=2)
    plt.plot(t0s,det.rate(t0s),'-', zorder=1)

plt.ylim(0)

plt.sca(ax[1])
plt.ylabel('z (t)')
for det in detectors:
    plt.plot(t0s, det.zs_c, '-', label='Counting analysis')

plt.sca(ax[2])
plt.ylabel('z (t)')
for det in detectors:
    #plt.plot(t0s, det.zs_c, ':', label='Counting analysis')
    plt.plot(t0s, det.zs_s, '-', label='Shape analysis')

plt.plot(t0s, zs_total, '-k', label='Shape analysis combined')
plt.legend()

plt.xlim(t0s[0],t0s[-1])
for a in ax:
    a.axvline(Tsn, color='k', lw=2)
    a.grid()

plt.xlabel('Time, s')
#plt.subplots_adjust(hspace=0.05)
plt.savefig('casa2.png')

