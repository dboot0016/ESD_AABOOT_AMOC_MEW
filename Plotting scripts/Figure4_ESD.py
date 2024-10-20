#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Earth System Dynamics
# Author script: A.A. Boot
# Contact person: Amber Boot (she/they; a.a.boot@uu.nl)

# Script for plotting Figure 4

#%% Import modules
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#%% Directory with data
Datadir='/Users/Boot0016/Documents/Stuff/Desktop/Cimatoribus_ESD_version/Data'

#%% Plotting variables
varx = 'E$_a$' # Variable on x-axis plots
vary = 'CO$_2$' # Variable on y-axis plots; For AMOC plots: 'D', for CO2 plots: 'CO$_2$'
unitx = '[Sv]' # Unit on x-axis
vers = '2' # What Es-coupling equation to use: '2' or '3'

save_fig = 'yes' # Whether to save figure: 'yes' to save
quality = 300 # Quality figure

LW = 4 # Linewidth of line
MS = 150 # Markersize
FS = 20 # Base fontsize in plots

#%% Model constants
D_dim = 1000 # Dimensionalize pycnocline depth
S_dim = 35 # Dimensionalize salinity

eta = 3e4 # Hydraulic constant
alpha = 2e-4 # Thermal expansion coefficient
beta = 8e-4 # Haline contraction coefficient
Tts = 10 # Temperature box ts
Tn = 5 # Temperature box n

tau = 0.1 # Wind stress
LxS = 3e7 # Zonal length Southern Ocean
rho0 = 1027.5 # Base density
fs = -1e-04 # Coriolis parameter
Agm = 1700 # Eddy diffusivity 
LxA = 1e7 # Zonal length Atlantic basin
Ly = 1e6 # Meridional extent frontal zone Southern Ocean

V0 = 3e17 # Total volume domain      
Vn = 3e15 # Volume box n       
Vs = 9e15 # Volume box s
A = 1e14 # Surface area box t 

epsilon = 1e-12 # Parameter for heavisede function
heaviside = lambda x: 0.5*(np.tanh(x/epsilon)+1) # Heaviside function

#%% Variables to load in data
a1 = np.repeat(0,13)
a2 = ([0,0,1,1,2,2,3,3,4,4,5,5,6])
a3 = np.tile([0,5],6)
a3 = np.append(a3,0)

#%% Create empty arrays for saddle node locations
L1 = np.zeros((13,))
L2 = np.zeros((13,))

#%% Select saddle node locations
for i in range(13):
    run1 = pd.read_csv(f'{Datadir}/'+'b.uncoupled_'+str(int(a1[i]))+str(a2[i])+str(a3[i])+'_branch_1',sep='\s+',header=None)   
    
    EA1 = run1[4]
    
    L1[i] = max(EA1)
    L2[i] = np.array(EA1)[-1]

#%% Es-CO2 couplings
Es = np.arange(0,0.61,0.05)

Es2 = np.arange(0,0.61,0.005)
CO2_2 = np.exp(((Es2+0.142)/0.092))

#%% Figure 4
new_tick_locations = np.array([0.25, 0.3,0.35,0.4,0.45,0.50])

def tick_function(X):
    V = np.exp(((X+0.142)/0.092))
    return ["%.0f" % z for z in V]

plt.figure(figsize=(8, 8))
fig, ax1 = plt.subplots(figsize=(8, 7.25))

plt.plot(Es,L1,'-.',color='tab:blue',linewidth=LW,label='On branch')
plt.plot(Es,L2,'--',color='tab:orange',linewidth=LW,label='Off branch')

plt.xlim([0.21,0.54])
plt.xlabel('E$_s$ [Sv]',fontsize=FS-2)
plt.ylabel('E$_a$ [Sv]',color='black',fontsize=FS-2)

plt.xticks([0.25,0.3,0.35,0.4,0.45,0.5],fontsize=FS-4)
plt.yticks([-0.15,0,0.15,0.3,0.45],color='black',fontsize=FS-4)    
plt.grid()
plt.legend(fontsize=FS-3,handlelength=3,loc=6)

plt.fill_between(Es, L1, 0.6, color='tab:blue', alpha = 0.2)
plt.fill_between(Es, L2, -0.25, color='tab:orange', alpha = 0.2)

plt.text(0.225,0.44,'Monostable off',fontsize = FS)
plt.text(0.225,-0.15,'Monostable on',fontsize = FS) 
plt.text(0.47,0.175,'MEW',fontsize = FS) 
plt.arrow(0.52,-0.04,0.0,0.45,color='black',width = 0.003,head_length = 0.03,shape='full')
plt.arrow(0.52,0.41,0.0,-0.44,color='black',width = 0.003,head_length = 0.03,shape='full')
plt.ylim([-0.2,0.5])

ax2 = ax1.twiny()
ax2.set_xlabel('CO$_2$ [ppm]', color='k',fontsize=FS-2)
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations),fontsize=FS-4)

plt.xlim([0.21,0.54])

if save_fig == 'yes':
    plt.savefig('figure_4.png', format='png', dpi=quality, bbox_inches = 'tight')
    
