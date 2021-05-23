#!/bin/env python 

import sn_stat as sn
import numpy as np
import pylab as plt

#define the rates
S = sn.rate(([0,1,10],[0,2,0]))*100
B = sn.rate(800)

Tsn = 0.
#total event rate vs time
R = B+S.shift(Tsn)

T0,T1 = -60,60

#create the dataset via sampler
s = sn.Sampler(R, time_window=[T0,T1])
ts = s.sample()

#detector configuration
det = sn.DetConfig(B=B, S=S)
print(det)

#assumed supernova start times
t0s = np.linspace(T0+det.time_window[0], T1-det.time_window[1],501)

sa = sn.ShapeAnalysis ([det])
ca = sn.CountingAnalysis(det)

#calculate TS vs time
ls_s = sa.l_val([ts],t0s)
ls_c = ca.l_val(ts,t0s)

#calculate z vs time
zs_s = sa([ts],t0s)
zs_c = ca(ts,t0s)

#plot everything
plt.style.use('seaborn-talk')

nplots = 3
fig, ax = plt.subplots(nplots,1, figsize=(12,4*nplots), sharex=True, constrained_layout=True)
plt.sca(ax[0])
plt.hist(ts,bins=np.arange(T0,T1,1), density=False)
plt.plot(t0s,R(t0s),'k-', zorder=1)
plt.ylim(0)
plt.ylabel('Event rate, 1/s')
plt.twinx()
plt.scatter(x=ts, y=np.random.uniform(size=len(ts)), marker='.', alpha=0.5, zorder=2, color='k', s=3)
plt.gca().get_yaxis().set_visible(False)

plt.sca(ax[1])
hc = plt.plot(t0s, ls_c, '-k', label='Counting analysis')[0]
plt.ylabel('N (t)')

plt.twinx()
hs = plt.plot(t0s, ls_s, '-b', label='Shape analysis')[0]
plt.ylabel('LLR (t)')
plt.legend(handles=[hc,hs], labels=['Counting analysis','Shape analysis'])

plt.sca(ax[2])
plt.plot(t0s, zs_c, '-k', label='Counting analysis')
plt.plot(t0s, zs_s, '-b', label='Shape analysis')
plt.ylabel('z (t)')
plt.legend()

plt.xlim(t0s[0],t0s[-1])
for a in ax:
    a.axvline(Tsn, color='r', lw=2)

plt.xlabel('Time, s')
#plt.subplots_adjust(hspace=0.05)
plt.savefig('llr.png')

