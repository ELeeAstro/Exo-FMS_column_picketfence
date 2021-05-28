import numpy as np
import matplotlib.pylab as plt
import pandas as pd

data1 = np.loadtxt('FMS_RC_ic.out')
Pic = data1[:,1]
Tic = data1[:,2]

data2 = np.loadtxt('FMS_RC_pp.out')
Ppp = data2[:,1]
Tpp = data2[:,2]
dT_rad = data2[:,3]
dT_conv = data2[:,4]
tauV = data2[:,5:8]
tauIR = data2[:,8:10]
kV = data2[:,10:13]
kIR = data2[:,13:15]


fig = plt.figure()

plt.plot(Tic,Pic/1e5,ls='dashed',lw=3,label='Initial Conditions',c='orange')
plt.plot(Tpp,Ppp/1e5,lw=3,label='Numerical Result',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('Temperature [K]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(dT_rad,Ppp/1e5,lw=3,label='dT_rad',c='red')
plt.plot(dT_conv,Ppp/1e5,lw=3,label='dT_conv',c='blue')
plt.ylabel('Pressure [bar]')
plt.xlabel('dT [K s$^{-1}$]')

plt.legend()
plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(tauV[:,0],Ppp/1e5,lw=3,label='tauV_1',c='midnightblue')
plt.plot(tauV[:,1],Ppp/1e5,lw=3,label='tauV_2',c='blue')
plt.plot(tauV[:,2],Ppp/1e5,lw=3,label='tauV_3',c='teal')
plt.plot(tauIR[:,0],Ppp/1e5,lw=3,label='tauIR_1',c='red')
plt.plot(tauIR[:,1],Ppp/1e5,lw=3,label='tauIR_2',c='darkorange')

plt.ylabel('Pressure [bar]')
plt.xlabel('Optical Depth')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(kV[:,0],Ppp/1e5,lw=3,label='kV_1',c='midnightblue')
plt.plot(kV[:,1],Ppp/1e5,lw=3,label='kV_2',c='blue')
plt.plot(kV[:,2],Ppp/1e5,lw=3,label='kV_3',c='teal')
plt.plot(kIR[:,0],Ppp/1e5,lw=3,label='kIR_1',c='red')
plt.plot(kIR[:,1],Ppp/1e5,lw=3,label='kIR_2',c='darkorange')

plt.ylabel('Pressure [bar]')
plt.xlabel('Opacity [m$^{2}$ kg$^{-1}$] ')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()


plt.show()
