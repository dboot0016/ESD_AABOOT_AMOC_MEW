#%% 

# Title paper: 'Potential effects of the marine carbon cycle on the multiple equilibria window of the Atlantic Meridional Overturning Circulation
# Authors: Boot A., von der Heydt, A.S., and Dijkstra, H.A.
# Submitted to Earth System Dynamics
# Author script: A.A. Boot
# Contact person: Amber Boot (she/they; a.a.boot@uu.nl)

# Script for plotting Figure S1

#%% Import modules
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#%% Directory with data
Datadir='/Users/daan/Desktop/Cimatoribus_ESD_version/Model_code/'

#%% Plotting variables
varx = 'E$_a$' # Variable on x-axis plots
vary = 'D' # Variable on y-axis plots; For AMOC plots: 'D', for CO2 plots: 'CO$_2$'
unitx = '[Sv]' # Unit on x-axis
vers = '2' # What Es-coupling equation to use: '2' or '3'

save_fig = 'yes' # Whether to save figure: 'yes' to save
quality = 300 # Quality figure

LW = 4 # Linewidth of line
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

#%% Function to load in variables from dataset
def ea_amoc_co2(run1,run2):
    EA1 = run1[4]
    St1 = run1[6]*S_dim
    Sts1 = run1[7]*S_dim
    Sn1 = run1[8]*S_dim
    Ss1 = run1[9]*S_dim
    D1 = run1[10]*D_dim
    CO21 = run1[23]*1e6
    
    Sd1 = (S_dim*V0 - St1*(A*D1)-Sts1*(LxA*Ly*D1)-Ss1*Vs-Sn1*Vn)/(V0-Vn-Vs-LxA*Ly*D1-A*D1)
    qn1 = 1e-6*(eta*alpha*(Tts-Tn)*(D1)**2+eta*beta*(Sn1-Sts1)*(D1)**2)
    AMOC1 = qn1*heaviside(qn1)
    
    EA2 = run2[4]
    St2 = run2[6]*S_dim
    Sts2 = run2[7]*S_dim
    Sn2 = run2[8]*S_dim
    Ss2 = run2[9]*S_dim
    D2 = run2[10]*D_dim
    CO22 = run2[23]*1e6
    
    Sd2 = (S_dim*V0 - St2*(A*D2)-Sts2*(LxA*Ly*D2)-Ss2*Vs-Sn2*Vn)/(V0-Vn-Vs-LxA*Ly*D2-A*D2)
    qn2 = 1e-6*(eta*alpha*(Tts-Tn)*(D2)**2+eta*beta*(Sn2-Sts2)*(D2)**2)
    AMOC2 = qn2*heaviside(qn2)
    
    return EA1,CO21,AMOC1,EA2,CO22,AMOC2

#%% Function for plotting results
def plot(x,y,run1,color):
    test = np.zeros((np.size(x),))
    for i in np.arange(0,np.size(x),1):
        if np.array(run1[1])[i]>0:
            test[i] = 1
        else:
           test[i] = 0
    
    test_size=np.size(np.where(abs(np.diff(test))==1))
    if test_size > 1:
        print('TEST IS LARGER THAN 1: CHECK STABILITY x1')
        print(test)
        
    plt.plot(x[:int(np.where(np.diff(test)==1)[0])+2],y[:int(np.where(np.diff(test)==1)[0])+2],linewidth=LW,color=color,linestyle = 'solid',label='_nolegend_')
    plt.plot(x[int(np.where(np.diff(test)==1)[0])+1:],y[int(np.where(np.diff(test)==1)[0])+1:],linewidth=LW,color=color,linestyle = 'dotted',label='_nolegend_')

#%% Load in data
run2_1 = pd.read_csv(f'{Datadir}/'+'b.es2_branch_1',sep='\s+',header=None)   
run1_1 = pd.read_csv(f'{Datadir}/'+'b.bio_branch_1',sep='\s+',header=None) 
run1_2 = pd.read_csv(f'{Datadir}/'+'b.bio_branch_2',sep='\s+',header=None)
run2_2 = pd.read_csv(f'{Datadir}/'+'b.es2_branch_2',sep='\s+',header=None)

for jj in range(5): 
    # Branch 1
    run3_1 = pd.read_csv(f'{Datadir}/'+'b.s'+str(jj+1)+'_branch_1',sep='\s+',header=None)   
    run4_1 = pd.read_csv(f'{Datadir}/'+'b.s'+str(jj+6)+'_branch_1',sep='\s+',header=None)    
    
    # Branch 2  
    run3_2 = pd.read_csv(f'{Datadir}/'+'b.s'+str(jj+1)+'_branch_2',sep='\s+',header=None)   
    run4_2 = pd.read_csv(f'{Datadir}/'+'b.s'+str(jj+6)+'_branch_2',sep='\s+',header=None)   
    
    #%% Variables
    EA1_1,CO21_1,AMOC1_1,EA1_2,CO21_2,AMOC1_2 = ea_amoc_co2(run1_1,run1_2)
    EA2_1,CO22_1,AMOC2_1,EA2_2,CO22_2,AMOC2_2 = ea_amoc_co2(run2_1,run2_2)
    EA3_1,CO23_1,AMOC3_1,EA3_2,CO23_2,AMOC3_2 = ea_amoc_co2(run3_1,run3_2)
    EA4_1,CO24_1,AMOC4_1,EA4_2,CO24_2,AMOC4_2 = ea_amoc_co2(run4_1,run4_2)
    
    #%% 
    ymin = 230 # minimum y-axis
    ymax = 350 # maximum y-axis
    ymax2 = 350 # maximum y-axis
    unity = '[ppm]' # unit y-axis
    vary2 = 'CO$_2$' # label y-axis
    
    #%% Fig S1b
    fig = plt.figure(figsize=(8, 6))

    label1='BIO'
    label2='BIO + E$_s$'
    label3='S'+(str(jj+1))
    label4='S'+str(jj+6)
    
    plot(EA1_1,CO21_1,run1_1,'black')
    plot(EA1_2,CO21_2,run1_2,'black')
    plot(EA2_1,CO22_1,run2_1,'tab:orange')
    plot(EA2_2,CO22_2,run2_2,'tab:orange')
    plot(EA3_1,CO23_1,run3_1,'tab:blue')
    plot(EA3_2,CO23_2,run3_2,'tab:blue')
    plot(EA4_1,CO24_1,run4_1,'tab:green')
    plot(EA4_2,CO24_2,run4_2,'tab:green')
    
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='black',label=label1)
    plt.plot([max(EA1_1),max(EA1_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='black',zorder=1)
    plt.plot([min(EA1_2),min(EA1_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='black',zorder=1)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:orange',label=label2)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:blue',label=label3)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:green',label=label4)
    
    plt.plot([min(EA1_2),min(EA1_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:orange',zorder=1)
    plt.plot([max(EA1_1),max(EA1_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:orange',zorder=1)
    
    plt.plot([min(EA2_2),min(EA2_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:blue',zorder=1)
    plt.plot([max(EA2_1),max(EA2_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:blue',zorder=1)
    
    plt.plot([min(EA3_2),min(EA3_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:blue',zorder=1)
    plt.plot([max(EA3_1),max(EA3_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:blue',zorder=1)
    
    plt.plot([min(EA4_2),min(EA4_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:green',zorder=1)
    plt.plot([max(EA4_1),max(EA4_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:green',zorder=1)
    
    plt.xlabel(str(varx)+' ' +str(unitx),fontsize=FS-2)
    plt.ylabel(str(vary2)+' ' +str(unity),fontsize=FS-2)
    
    plt.xticks(fontsize=FS-4)
    plt.yticks(fontsize=FS-4)    
    plt.grid()

    plt.xlim([-0.22,0.5])
    plt.ylim([ymin,ymax2])
    
    if save_fig == 'yes':
        plt.savefig('figureS1_CO2_'+str(jj)+'.png', format='png', dpi=quality, bbox_inches = 'tight')

    #%%
    ymin = -0.75 # minimum y-axis
    ymax = 23 # maximum y-axis
    ymax2 = 23 # maximum y-axis
    unity = '[Sv]' # unit y-axis
    vary2 = 'AMOC' # label y-axis
 
    #%% Figure S1a   
    fig = plt.figure(figsize=(8, 6))

    label1='BIO'
    label2='BIO + E$_s$'
    label3='S'+(str(jj+1))
    label4='S'+str(jj+6)
    
    plot(EA1_1,AMOC1_1,run1_1,'tab:blue')
    plot(EA1_2,AMOC1_2,run1_2,'tab:blue')
    plot(EA2_1,AMOC2_1,run2_1,'black')
    plot(EA2_2,AMOC2_2,run2_2,'black')
    plot(EA3_1,AMOC3_1,run3_1,'tab:orange')
    plot(EA3_2,AMOC3_2,run3_2,'tab:orange')
    plot(EA4_1,AMOC4_1,run4_1,'tab:green')
    plot(EA4_2,AMOC4_2,run4_2,'tab:green')
        
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='black',label=label1)
    plt.plot([max(EA1_1),max(EA1_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='black',zorder=1)
    plt.plot([min(EA1_2),min(EA1_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='black',zorder=1)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:orange',label=label2)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:blue',label=label3)
    plt.plot([-1e5,-2e5],[np.mean(CO21_1),np.mean(CO24_1)],linewidth=LW,color='tab:green',label=label4)
    
    # Saddle node locations
    plt.plot([min(EA2_2),min(EA2_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:orange',zorder=1)
    plt.plot([max(EA2_1),max(EA2_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:orange',zorder=1)
    plt.plot([min(EA3_2),min(EA3_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:blue',zorder=1)
    plt.plot([max(EA3_1),max(EA3_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:blue',zorder=1)
    plt.plot([min(EA4_2),min(EA4_2)],[-10,ymax],linewidth=LW-2,linestyle = 'dashed',color='tab:green',zorder=1)
    plt.plot([max(EA4_1),max(EA4_1)],[-10,ymax],linewidth=LW-2,linestyle = 'dashdot',color='tab:green',zorder=1)
    
    plt.xlabel(str(varx)+' ' +str(unitx),fontsize=FS-2)
    plt.ylabel(str(vary2)+' ' +str(unity),fontsize=FS-2)
    
    plt.xticks(fontsize=FS-4)
    plt.yticks(fontsize=FS-4)    
    plt.grid()
    plt.legend(fontsize=FS-6)
    
    plt.xlim([-0.22,0.5])
    plt.ylim([ymin,ymax2])
    
    if save_fig == 'yes':
        plt.savefig('figureS1_AMOC_'+str(jj)+'.png', format='png', dpi=quality, bbox_inches = 'tight')
