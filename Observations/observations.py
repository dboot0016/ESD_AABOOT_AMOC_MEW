#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import xesmf as xe
import cmocean as cm
import regionmask
from cmip6_preprocessing.regionmask import merged_mask

from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import scipy

#%% Load data
datadir='/Users/daan/Downloads/ORD49465/' 

load_var1a = xr.open_dataset(f'{datadir}/EMPmm19870701000000213SCPOS01UD.nc' )
    
yrs = np.arange(1987,2015,1)
mon = (['01','02','03','04','05','06','07','08','09','10','11','12'])

yrs = np.repeat(yrs,12)
yrs = yrs[7:]
mon = np.tile(mon,np.size(yrs/12))
mon = mon[7:]

#%%
for i in range(len(yrs)):
    load_var1b = xr.open_dataset(f'{datadir}/EMPmm'+str(yrs[i])+mon[i]+'01000000213SCPOS01UD.nc' )
    load_var1a = xr.concat([load_var1a,load_var1b],dim='time')

#%% Weighted temporal mean function
def weighted_temporal_mean(ds, loc, dt):
    month_length = ds.time.dt.days_in_month
    wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
    np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)
    
    obs = ds[loc]
    cond = obs.isnull()
    ones = xr.where(cond, 0.0, 1.0)

    obs_sum = (obs * wgts).resample(time=str(dt)+"AS").sum(dim="time")
    ones_out = (ones * wgts).resample(time=str(dt)+"AS").sum(dim="time")

    return obs_sum / ones_out

#%%
var = weighted_temporal_mean(load_var1a,'budg',1)*1e-3/86400

#%%
lat=var.lat
lon=var.lon

RE = 6.371e6  # [m] Earth radius
                             # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360
dy = -2*np.pi*RE*(lat[1]-lat[0]).values/360       

#%%
var1 = var
var1[:,-50:,:60] = 'nan'
var1[:,-78:-50,:50] = 'nan'
var1[:,68:-78,:30] = 'nan'
var1[:,:50,198:] = 'nan'
#%%
VAR_lon=(var1*dx).sum(['lon'])            # Averaged over lon
VAR_lat=(VAR_lon[:,:]*dy).sum(['lat'])                   # Averaged over lat
VAR1=VAR_lat
