#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Earth System Dynamics
# Author script: A.A. Boot
# Contact person: Amber Boot (she/they; a.a.boot@uu.nl)

# Script for plotting Figure 5

#%% Import modules
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#%% Directory with data
Datadir='/Users/daan/Desktop/Cimatoribus_ESD_version/Model_code/'

#%% Plotting variables
varx = 'E$_a$' # Variable on x-axis plots
vary = 'CO$_2$' # Variable on y-axis plots; For AMOC plots: 'D', for CO2 plots: 'CO$_2$'
unitx = '[Sv]' # Unit on x-axis
vers = '2' # What Es-coupling equation to use: '2' or '3'

save_fig = 'on' # Whether to save figure: 'yes' to save
quality = 300 # Quality figure

LW = 3.2 # Linewidth of line
MS = 60 # Markersize
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

#%% Set up plotting domain
y1 = 0.995 + 0.001
y2 = 0.995 + 0.003
y3 = 0.995 + 0.005
y4 = 0.995 + 0.007
y5 = 0.995 + 0.009

y6 = 0.995 + 0.011
y7 = 0.995 + 0.013
y8 = 0.995 + 0.015
y9 = 0.995 + 0.017
y10 = 0.995 + 0.019

y11 = 0.995 + 0.021
y12 = 0.995 + 0.023
y13 = 0.995 + 0.025
y14 = 0.995 + 0.027
y15 = 0.995 + 0.029

y16 = 0.995 + 0.031
y17 = 0.995 + 0.033
y18 = 0.995 + 0.035
y19 = 0.995 + 0.037
y20 = 0.995 + 0.039

y21 = 0.995 + 0.041
y22 = 0.995 + 0.043
y23 = 0.995 + 0.045
y24 = 0.995 + 0.047
y25 = 0.995 + 0.049

#%% CS_HI
run5_1 = pd.read_csv(f'{Datadir}/'+'b.temp22_m4000_branch_1',sep='\s+',header=None)  
run4_1 = pd.read_csv(f'{Datadir}/'+'b.temp22_branch_1',sep='\s+',header=None)  
run3_1 = pd.read_csv(f'{Datadir}/'+'b.temp22_4000_branch_1',sep='\s+',header=None)  
run2_1 = pd.read_csv(f'{Datadir}/'+'b.temp22_10000_branch_1',sep='\s+',header=None)  
run1_1 = pd.read_csv(f'{Datadir}/'+'b.temp22_20000_branch_1',sep='\s+',header=None)

run5_2 = pd.read_csv(f'{Datadir}/'+'b.temp22_m4000_branch_2',sep='\s+',header=None)  
run4_2 = pd.read_csv(f'{Datadir}/'+'b.temp22_branch_2',sep='\s+',header=None)  
run3_2 = pd.read_csv(f'{Datadir}/'+'b.temp22_4000_branch_2',sep='\s+',header=None)  
run2_2 = pd.read_csv(f'{Datadir}/'+'b.temp22_10000_branch_2',sep='\s+',header=None)  
run1_2 = pd.read_csv(f'{Datadir}/'+'b.temp22_20000_branch_2',sep='\s+',header=None)  

#%% CS_LO
run10_1 = pd.read_csv(f'{Datadir}/'+'b.temp12_m4000_branch_1',sep='\s+',header=None)  
run9_1 = pd.read_csv(f'{Datadir}/'+'b.temp12_branch_1',sep='\s+',header=None)  
run8_1 = pd.read_csv(f'{Datadir}/'+'b.temp12_4000_branch_1',sep='\s+',header=None)  
run7_1 = pd.read_csv(f'{Datadir}/'+'b.temp12_10000_branch_1',sep='\s+',header=None)  
run6_1 = pd.read_csv(f'{Datadir}/'+'b.temp12_20000_branch_1',sep='\s+',header=None)

run10_2 = pd.read_csv(f'{Datadir}/'+'b.temp12_m4000_branch_2',sep='\s+',header=None)  
run9_2 = pd.read_csv(f'{Datadir}/'+'b.temp12_branch_2',sep='\s+',header=None)  
run8_2 = pd.read_csv(f'{Datadir}/'+'b.temp12_4000_branch_2',sep='\s+',header=None)  
run7_2 = pd.read_csv(f'{Datadir}/'+'b.temp12_10000_branch_2',sep='\s+',header=None)  
run6_2 = pd.read_csv(f'{Datadir}/'+'b.temp12_20000_branch_2',sep='\s+',header=None)  

#%% FCA
run15_1 = pd.read_csv(f'{Datadir}/'+'b.fca2_m4000_branch_1',sep='\s+',header=None)  
run14_1 = pd.read_csv(f'{Datadir}/'+'b.fca2_branch_1',sep='\s+',header=None)  
run13_1 = pd.read_csv(f'{Datadir}/'+'b.fca2_4000_branch_1',sep='\s+',header=None)  
run12_1 = pd.read_csv(f'{Datadir}/'+'b.fca2_10000_branch_1',sep='\s+',header=None)  
run11_1 = pd.read_csv(f'{Datadir}/'+'b.fca2_20000_branch_1',sep='\s+',header=None)

run15_2 = pd.read_csv(f'{Datadir}/'+'b.fca2_m4000_branch_2',sep='\s+',header=None)  
run14_2 = pd.read_csv(f'{Datadir}/'+'b.fca2_branch_2',sep='\s+',header=None)  
run13_2 = pd.read_csv(f'{Datadir}/'+'b.fca2_4000_branch_2',sep='\s+',header=None)  
run12_2 = pd.read_csv(f'{Datadir}/'+'b.fca2_10000_branch_2',sep='\s+',header=None)  
run11_2 = pd.read_csv(f'{Datadir}/'+'b.fca2_20000_branch_2',sep='\s+',header=None)  

#%% ES
run20_1 = pd.read_csv(f'{Datadir}/'+'b.es2_m4000_branch_1',sep='\s+',header=None)  
run19_1 = pd.read_csv(f'{Datadir}/'+'b.es2_branch_1',sep='\s+',header=None)  
run18_1 = pd.read_csv(f'{Datadir}/'+'b.es2_4000_branch_1',sep='\s+',header=None)  
run17_1 = pd.read_csv(f'{Datadir}/'+'b.es2_10000_branch_1',sep='\s+',header=None)  
run16_1 = pd.read_csv(f'{Datadir}/'+'b.es2_20000_branch_1',sep='\s+',header=None)

run20_2 = pd.read_csv(f'{Datadir}/'+'b.es2_m4000_branch_2',sep='\s+',header=None)  
run19_2 = pd.read_csv(f'{Datadir}/'+'b.es2_branch_2',sep='\s+',header=None)  
run18_2 = pd.read_csv(f'{Datadir}/'+'b.es2_4000_branch_2',sep='\s+',header=None)  
run17_2 = pd.read_csv(f'{Datadir}/'+'b.es2_10000_branch_2',sep='\s+',header=None)  
run16_2 = pd.read_csv(f'{Datadir}/'+'b.es2_20000_branch_2',sep='\s+',header=None)

#%% BIO
run25_1 = pd.read_csv(f'{Datadir}/'+'b.bio_m4000_branch_1',sep='\s+',header=None)  
run24_1 = pd.read_csv(f'{Datadir}/'+'b.bio_branch_1',sep='\s+',header=None)  
run23_1 = pd.read_csv(f'{Datadir}/'+'b.bio_4000_branch_1',sep='\s+',header=None)  
run22_1 = pd.read_csv(f'{Datadir}/'+'b.bio_10000_branch_1',sep='\s+',header=None)  
run21_1 = pd.read_csv(f'{Datadir}/'+'b.bio_20000_branch_1',sep='\s+',header=None)

run25_2 = pd.read_csv(f'{Datadir}/'+'b.bio_m4000_branch_2',sep='\s+',header=None)  
run24_2 = pd.read_csv(f'{Datadir}/'+'b.bio_branch_2',sep='\s+',header=None)  
run23_2 = pd.read_csv(f'{Datadir}/'+'b.bio_4000_branch_2',sep='\s+',header=None)  
run22_2 = pd.read_csv(f'{Datadir}/'+'b.bio_10000_branch_2',sep='\s+',header=None)  
run21_2 = pd.read_csv(f'{Datadir}/'+'b.bio_20000_branch_2',sep='\s+',header=None)  

#%% Select Ea data
def select_plot_data(run1,run2):
    EA1 = run1[4]
    EA2 = run2[4]
    
    CO21 = run1[23]*1e6
    CO22 = run2[23]*1e6
    
    X1_1 = max(EA1)
    X1_2 = min(EA2)
    
    X2_1 = CO21[np.where(max(EA1)==EA1)[0][0]]
    X2_2 = CO22[np.where(min(EA2)==EA2)[0][0]]
    
    ln_a = -0.142
    ln_b = 0.097
    ln_c = 0
    ln_d = 1

    X3_1 = ln_a+(ln_b*np.log((X2_1+ln_c)/ln_d))
    X3_2 = ln_a+(ln_b*np.log((X2_2+ln_c)/ln_d))
    
    return X1_1,X1_2,X2_1,X2_2,X3_1,X3_2

#%% Make plot domain
def plot_base():
    Y=np.zeros(100)+1
    xline=np.linspace(-1500,1500,np.size(Y))
    
    X_out = np.linspace(10,11,100)
    Y_out = np.linspace(10,11,100)
    
    plt.plot(xline,Y+0.005,'k',linewidth=0.75)
    plt.plot(xline,Y+0.015,'k',linewidth=0.75)
    plt.plot(xline,Y+0.025,'k',linewidth=0.75)
    plt.plot(xline,Y+0.035,'k',linewidth=0.75)
    plt.plot(xline,Y+0.045,'k',linewidth=0.75)

#%% Plot data
def plot_data(x1,y1,x2,color):
    plt.scatter(x1,y1,s=MS,color=color,label='_nolegend_')
    plt.scatter(x2,y1,s=MS,marker='s',color=color,label='_nolegend_')
    plt.plot([x1,x2],[y1,y1],linewidth=LW,linestyle='dashed',color=color,label='_nolegend_')
    
#%% Select data
X1_11,X1_21,X1_12,X1_22,X1_13,X1_23 = select_plot_data(run1_1,run1_2)
X2_11,X2_21,X2_12,X2_22,X2_13,X2_23 = select_plot_data(run2_1,run2_2)
X3_11,X3_21,X3_12,X3_22,X3_13,X3_23 = select_plot_data(run3_1,run3_2)
X4_11,X4_21,X4_12,X4_22,X4_13,X4_23 = select_plot_data(run4_1,run4_2)
X5_11,X5_21,X5_12,X5_22,X5_13,X5_23 = select_plot_data(run5_1,run5_2)

X6_11,X6_21,X6_12,X6_22,X6_13,X6_23 = select_plot_data(run6_1,run6_2)
X7_11,X7_21,X7_12,X7_22,X7_13,X7_23 = select_plot_data(run7_1,run7_2)
X8_11,X8_21,X8_12,X8_22,X8_13,X8_23 = select_plot_data(run8_1,run8_2)
X9_11,X9_21,X9_12,X9_22,X9_13,X9_23 = select_plot_data(run9_1,run9_2)
X10_11,X10_21,X10_12,X10_22,X10_13,X10_23 = select_plot_data(run10_1,run10_2)

X11_11,X11_21,X11_12,X11_22,X11_13,X11_23 = select_plot_data(run11_1,run11_2)
X12_11,X12_21,X12_12,X12_22,X12_13,X12_23 = select_plot_data(run12_1,run12_2)
X13_11,X13_21,X13_12,X13_22,X13_13,X13_23 = select_plot_data(run13_1,run13_2)
X14_11,X14_21,X14_12,X14_22,X14_13,X14_23 = select_plot_data(run14_1,run14_2)
X15_11,X15_21,X15_12,X15_22,X15_13,X15_23 = select_plot_data(run15_1,run15_2)

X16_11,X16_21,X16_12,X16_22,X16_13,X16_23 = select_plot_data(run16_1,run16_2)
X17_11,X17_21,X17_12,X17_22,X17_13,X17_23 = select_plot_data(run17_1,run17_2)
X18_11,X18_21,X18_12,X18_22,X18_13,X18_23 = select_plot_data(run18_1,run18_2)
X19_11,X19_21,X19_12,X19_22,X19_13,X19_23 = select_plot_data(run19_1,run19_2)
X20_11,X20_21,X20_12,X20_22,X20_13,X20_23 = select_plot_data(run20_1,run20_2)

X21_11,X21_21,X21_12,X21_22,X21_13,X21_23 = select_plot_data(run21_1,run21_2)
X22_11,X22_21,X22_12,X22_22,X22_13,X22_23 = select_plot_data(run22_1,run22_2)
X23_11,X23_21,X23_12,X23_22,X23_13,X23_23 = select_plot_data(run23_1,run23_2)
X24_11,X24_21,X24_12,X24_22,X24_13,X24_23 = select_plot_data(run24_1,run24_2)
X25_11,X25_21,X25_12,X25_22,X25_13,X25_23 = select_plot_data(run25_1,run25_2)

#%% Figure 5a
fig, ax = plt.subplots(figsize=(8,6))
plot_base()

plot_data(X1_11,y1,X1_21,'tab:blue')
plot_data(X2_11,y2,X2_21,'tab:orange')
plot_data(X3_11,y3,X3_21,'tab:green')
plot_data(X4_11,y4,X4_21,'black')
plot_data(X5_11,y5,X5_21,'tab:purple')

plot_data(X6_11,y6,X6_21,'tab:blue')
plot_data(X7_11,y7,X7_21,'tab:orange')
plot_data(X8_11,y8,X8_21,'tab:green')
plot_data(X9_11,y9,X9_21,'black')
plot_data(X10_11,y10,X10_21,'tab:purple')

plot_data(X11_11,y11,X11_21,'tab:blue')
plot_data(X12_11,y12,X12_21,'tab:orange')
plot_data(X13_11,y13,X13_21,'tab:green')
plot_data(X14_11,y14,X14_21,'black')
plot_data(X15_11,y15,X15_21,'tab:purple')

plot_data(X16_11,y16,X16_21,'tab:blue')
plot_data(X17_11,y17,X17_21,'tab:orange')
plot_data(X18_11,y18,X18_21,'tab:green')
plot_data(X19_11,y19,X19_21,'black')
plot_data(X20_11,y20,X20_21,'tab:purple')

plot_data(X21_11,y21,X21_21,'tab:blue')
plot_data(X22_11,y22,X22_21,'tab:orange')
plot_data(X23_11,y23,X23_21,'tab:green')
plot_data(X24_11,y24,X24_21,'black')
plot_data(X25_11,y25,X25_21,'tab:purple')

#plt.title('Bistability window',fontsize=FS)
plt.xticks(fontsize=FS-4)
plt.xlim([-0.12, 0.48])
plt.ylim([0.995, 1.045])
plt.yticks([1,1.01,1.02,1.03,1.04],['CS$_{HI}$','CS$_{LO}$','FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('E$_a$ [Sv]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")

if save_fig == 'on':
    plt.savefig('Figure4_a.png', format='png', dpi=quality, bbox_inches = 'tight')
    
#%% Figure 5b
fig, ax = plt.subplots(figsize=(8,6))
plot_base()

plot_data(X1_12,y1,X1_22,'tab:blue')
plot_data(X2_12,y2,X2_22,'tab:orange')
plot_data(X3_12,y3,X3_22,'tab:green')
plot_data(X4_12,y4,X4_22,'black')
plot_data(X5_12,y5,X5_22,'tab:purple')

plot_data(X6_12,y6,X6_22,'tab:blue')
plot_data(X7_12,y7,X7_22,'tab:orange')
plot_data(X8_12,y8,X8_22,'tab:green')
plot_data(X9_12,y9,X9_22,'black')
plot_data(X10_12,y10,X10_22,'tab:purple')

plot_data(X11_12,y11,X11_22,'tab:blue')
plot_data(X12_12,y12,X12_22,'tab:orange')
plot_data(X13_12,y13,X13_22,'tab:green')
plot_data(X14_12,y14,X14_22,'black')
plot_data(X15_12,y15,X15_22,'tab:purple')

plot_data(X16_12,y16,X16_22,'tab:blue')
plot_data(X17_12,y17,X17_22,'tab:orange')
plot_data(X18_12,y18,X18_22,'tab:green')
plot_data(X19_12,y19,X19_22,'black')
plot_data(X20_12,y20,X20_22,'tab:purple')

plot_data(X21_12,y21,X21_22,'tab:blue')
plot_data(X22_12,y22,X22_22,'tab:orange')
plot_data(X23_12,y23,X23_22,'tab:green')
plot_data(X24_12,y24,X24_22,'black')
plot_data(X25_12,y25,X25_22,'tab:purple')

#plt.title('Bistability window',fontsize=FS)
plt.xticks(fontsize=FS-4)
plt.xlim([120, 800])
plt.ylim([0.995, 1.045])
plt.yticks([1,1.01,1.02,1.03,1.04],['CS$_{HI}$','CS$_{LO}$','FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('CO$_2$ [ppm]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")

X_out = np.linspace(10,11,100)
Y_out = np.linspace(10,11,100)
plt.plot(X_out,Y_out,linewidth=LW,color='tab:purple',label = '-4000 PgC')
plt.plot(X_out,Y_out,linewidth=LW,color='black',label = '0 PgC')
plt.plot(X_out,Y_out,linewidth=LW,color='tab:green',label = '+4000 PgC')
plt.plot(X_out,Y_out,linewidth=LW,color='tab:orange',label = '+10000 PgC')
plt.plot(X_out,Y_out,linewidth=LW,color='tab:blue',label = '+20000 PgC')
plt.legend(fontsize=FS-6,loc = 'upper right')

if save_fig == 'on':
    plt.savefig('Figure4_b.png', format='png', dpi=quality, bbox_inches = 'tight')

#%% Figure 5c
fig, ax = plt.subplots(figsize=(8,6))
plot_base()

plot_data(X1_13,y1,X1_23,'tab:blue')
plot_data(X2_13,y2,X2_23,'tab:orange')
plot_data(X3_13,y3,X3_23,'tab:green')
plot_data(X4_13,y4,X4_23,'black')
plot_data(X5_13,y5,X5_23,'tab:purple')

plot_data(X6_13,y6,X6_23,'tab:blue')
plot_data(X7_13,y7,X7_23,'tab:orange')
plot_data(X8_13,y8,X8_23,'tab:green')
plot_data(X9_13,y9,X9_23,'black')
plot_data(X10_13,y10,X10_23,'tab:purple')

plot_data(X11_13,y11,X11_23,'tab:blue')
plot_data(X12_13,y12,X12_23,'tab:orange')
plot_data(X13_13,y13,X13_23,'tab:green')
plot_data(X14_13,y14,X14_23,'black')
plot_data(X15_13,y15,X15_23,'tab:purple')

plot_data(X16_13,y16,X16_23,'tab:blue')
plot_data(X17_13,y17,X17_23,'tab:orange')
plot_data(X18_13,y18,X18_23,'tab:green')
plot_data(X19_13,y19,X19_23,'black')
plot_data(X20_13,y20,X20_23,'tab:purple')

plot_data(0.43,y21,0.43,'tab:blue')
plot_data(0.43,y22,0.43,'tab:orange')
plot_data(0.43,y23,0.43,'tab:green')
plot_data(0.43,y24,0.43,'black')
plot_data(0.43,y25,0.43,'tab:purple')

#plt.title('Bistability window',fontsize=FS)
plt.xticks(fontsize=FS-4)
plt.xlim([0.33, 0.52])
plt.xticks([0.35,0.38,0.41,0.44,0.47,0.50])
plt.ylim([0.995, 1.045])
plt.yticks([1,1.01,1.02,1.03,1.04],['CS$_{HI}$','CS$_{LO}$','FCA','E$_s$','BIO'],fontsize=FS-2)
plt.xlabel('E$_s$ [Sv]',fontsize=FS-2)
ax.xaxis.grid(True, which="minor", ls="--")
ax.xaxis.grid(True, which="major", ls="-")

plt.scatter(X_out,Y_out,s=MS,color='black',label = 'On')
plt.scatter(X_out,Y_out,s=MS,marker='s',color='black',label = 'Off')
plt.legend(fontsize=FS-6)

if save_fig == 'on':
    plt.savefig('Figure4_c.png', format='png', dpi=quality, bbox_inches = 'tight')