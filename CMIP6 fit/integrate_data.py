#%% 


# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A.A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Earth System Dynamics
# Author script: A.A. Boot
# Contact person: Amber Boot (she/they; a.a.boot@uu.nl)

# Script for integrating and resaving data
# Data necessary to run script can be downloaded from https://esgf-node.llnl.gov/search/cmip6/
# Data needs to be regridded using the concat_regrid_save.py script

#%% Import modules
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import xesmf as xe
import pandas as pd

#%% Data dir
datadir='/Users/daan/Desktop/ESD_MEW_Repository/' 
var1 = 'wfo'

#%% Area data
data1='/Users/daan/CESM2_data'   # Location of dataset(s) 
load_var1 = xr.open_dataset(f'{data1}/area_gn.nc')
area_gn=load_var1['areacello'][:,:].compute().squeeze()

load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
area_gr=load_var1['areacello'][:,:].compute().squeeze()

lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

load_var2 = xr.open_dataset(f'{datadir}/mask_atlantic_th.nc') # Mask for Atlantic thermocline region
mask_ATL=load_var2[var1].compute().squeeze()

load_var2 = xr.open_dataset(f'{datadir}/mask_pacific_th.nc') # Indo-Pacific mask
mask_PAC=load_var2[var1].compute().squeeze()
  
#%% Load in data
filenames_list = np.array(pd.read_csv(f'{datadir}/filenames_file.txt',sep='\s+',header=None)) # List with model names and conversion factos

def load_and_integrate(i):
    filename = 'wfo_'+str(filenames_list[i,0])+'_'+str(filenames_list[i,1])+'.nc'

    load_var1 = xr.open_dataset(f'{datadir}/'+filename)
    VAR1=load_var1[var1][:150*12,:,:].compute().squeeze()*filenames_list[i,2]

    V_ATL = (VAR1*mask_ATL*area_gr).sum(['lon','lat'])
    V_PAC = (VAR1*mask_PAC*area_gr).sum(['lon','lat']) 
    
    fig = plt.figure()
    plt.plot(V_ATL/1027/1e6,label = 'Atlantic')
    plt.plot(V_PAC/1027/1e6,label = 'Pacific')
    plt.title(str(filenames_list[i,0]))
    plt.legend()
    
    return V_ATL,V_PAC

#%% Load and integrate first data
atl1,pac1 = load_and_integrate(0)  
atl2,pac2 = load_and_integrate(1)

#%% Reshape data
atl = np.reshape(np.concatenate((atl1,atl2)),(1800,2),order='F')
pac = np.reshape(np.concatenate((pac1,pac2)),(1800,2),order='F') 

#%% Load in all data
for i in range(filenames_list.shape[0]-2):
    atl1,pac1 = load_and_integrate(i+2)
    atl = np.append(atl,np.array(atl1).reshape((1800,1)),axis=1)
    pac = np.append(pac,np.array(pac1).reshape((1800,1)),axis=1)
    print(i+1)                             
               
#%% Units to Sv                       
ATL = atl/1027/1e6 # Units in Sv
PAC = pac/1027/1e6 # Units in Sv

#%% Function for making a data array
def MakeDataArray(data, time, model, dims):
        
        array = xr.DataArray(data = data,
                                 dims = ["time", "model"],
                                 coords = dict(time = time,
                                               model = model));   
        return array
    
#%% Save data
datadir='/Users/daan/Desktop/ESD_MEW_Repository/' 
filename1 = 'ATL_ES_FW.nc'
filename2 = 'PAC_EP_FW.nc'

read_var1 = MakeDataArray(data = ATL, time=atl1['time'], model=filenames_list[:,0], dims = "1d")
read_var2 = MakeDataArray(data = PAC, time=atl1['time'], model=filenames_list[:,0], dims = "1d")

data1 = read_var1.to_dataset(name = 'wfo');
data1.to_netcdf(datadir+filename1)

data2 = read_var2.to_dataset(name = 'wfo');
data2.to_netcdf(datadir+filename2)



