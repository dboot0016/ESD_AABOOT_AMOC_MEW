#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Earth System Dynamics
# Author script: A.A. Boot
# Contact person: Amber Boot (she/they; a.a.boot@uu.nl)

# Script for plotting Figure 2

# Data is prepared using the following scripts:
    # - concat_regrid_save.py
    # - integrate_data.py

#%% Import modules
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import xesmf as xe
import pandas as pd

#%% Data dir
datadir='/Users/daan/Desktop/ESD_MEW_Repository/' 
var1 = 'wfo'

#%% Load data
load_var1 = xr.open_dataset(f'{datadir}/ATL_ES_FW.nc')
ATL=load_var1[var1].compute().squeeze()

load_var1 = xr.open_dataset(f'{datadir}/PAC_EP_FW.nc')
PAC=load_var1[var1].compute().squeeze()

#%% Create 5-yr running mean
tm = 60
Lv11 = pd.DataFrame(ATL)
rolling_windows = Lv11.rolling(tm)
ATL1=np.array(np.squeeze(rolling_windows.mean()))

Lv11 = pd.DataFrame(PAC)
rolling_windows = Lv11.rolling(tm)
PAC1=np.array(np.squeeze(rolling_windows.mean()))

#%% Load in CO2 dataset 
load_var1 = xr.open_dataset(f'{datadir}/co2mass_1pct_co2_CESM2_1.nc')
VAR1=load_var1['co2mass'].compute().squeeze()
CO2_1pct = VAR1/5.14e18*1.0e6 * 28.966 / 44.0

#%% Define fitting function
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

def func_lin(x,a,b):
    return a*x+b

def func_log(x,a,b):
    return a+b*np.log(x)

def fit(data):
    popt_1, pcov_1 = curve_fit(func_lin, np.array(CO2_1pct)[60:], np.array(data)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))
    popt_2, pcov_2 = curve_fit(func_log, np.array(CO2_1pct)[60:], np.array(data)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))

    fit1 = popt_1[0]*CO2_1pct+popt_1[1]
    fit2 = popt_2[0]+popt_2[1]*np.log(CO2_1pct)
    
    return fit1, fit2

#%% Create linear and log fit
data1=np.mean(ATL1,axis=1)/2
popt_1, pcov_1 = curve_fit(func_lin, np.array(CO2_1pct)[60:], np.array(data1)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))
popt_2, pcov_2 = curve_fit(func_log, np.array(CO2_1pct)[60:], np.array(data1)[60:],bounds=([-5e6,-5e6], [5e6, 5e6]))

fit1 = popt_1[0]*CO2_1pct+popt_1[1]
fit2 = popt_2[0]+popt_2[1]*np.log(CO2_1pct)

r21=r2_score(np.array(data1)[60:], fit1[60:])
r22=r2_score(np.array(data1)[60:], fit2[60:])

print('Linear fit:')
print(str(np.round(popt_1[1],3))+' + ' + str(np.round(popt_1[0],6)) +' x CO2')

print('Log fit:')
print(str(np.round(popt_2[0],3))+' + ' + str(np.round(popt_2[1],3)) +' x ln(CO2)')

#%% Figure 2c
LW = 4
MS = 50
FS = 20
save_fig = 'no'
quality = 300

fig = plt.figure(figsize=(7, 5)) 
plt.plot(CO2_1pct,data1,linewidth=LW,zorder=4,color='black',label='MM')
[fit1,fit2]=fit(data1)
#plt.plot(CO2_1pct,fit1,color='tab:orange',linewidth=LW,label='Linear fit; R$^2$ = ' +str(np.round(r21,2)))
plt.plot(CO2_1pct,fit2,color='tab:orange',linewidth=LW,label='Logarithmic fit; R$^2$ = ' +str(np.round(r22,2)))

plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
plt.title('Atlantic (multimodel mean)',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6)
plt.grid()

if save_fig == 'yes':
    plt.savefig('Figure2_c.png', format='png', dpi=quality, bbox_inches = 'tight')
 
#%% Figure 2v
fig = plt.figure(figsize=(7, 5))
for i in range(ATL1.shape[1]):
    plt.scatter(CO2_1pct,np.array(ATL1/2)[:,i],s=5)

plt.plot(CO2_1pct,data1,linewidth=LW,color='black',label='Multimodel mean')
plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_s$ [Sv]',fontsize=FS-2)
plt.title('Atlantic (CMIP6 ensemble)',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6,ncol=1,loc=2)
plt.grid()

if save_fig == 'yes':
    plt.savefig('Figure2_b.png', format='png', dpi=quality, bbox_inches = 'tight')
    
#%% Figure 2a
fig = plt.figure(figsize=(7, 5))
for i in range(PAC1.shape[1]):
    plt.scatter(CO2_1pct,np.array(PAC1)[:,i],s=5)

plt.plot(CO2_1pct,np.mean(PAC1,axis=1),linewidth=LW,color='black',label='Multimodel mean')
plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
plt.ylabel('E$_p$ [Sv]',fontsize=FS-2)
plt.title('Pacific (CMIP6 ensemble)',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(fontsize=FS-6,ncol=1,loc=2)
plt.grid()

if save_fig == 'yes':
    plt.savefig('Figure2_a.png', format='png', dpi=quality, bbox_inches = 'tight')

#%%
print('Mean Indo-Pacific flux: ' + str(np.nanmean(PAC1)) + ' Sv')

