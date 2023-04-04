# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 19:13:22 2022

@author: MOHANRAHAN - group 12 - 5627796
"""


import numpy as np
import pandas as pd
import scipy.integrate as spi
import matplotlib.pyplot as plt
import scipy.sparse as sp
import seaborn as sns
import MyTicToc as mt
import sys




sns.set()

####################################

'''Testing Scenarios can be set using the 1-3 index of the following list'''

run_scenario = 0

scenario = ['scenario 0',           #base case no root water uptake
            'Scenario 1',           #clay loam with volume to 30cm
            'Scenario 2',           #clay loam with Lrv to 2m 
            'Scenario 3']           #describe

####################################



meteo    = pd.read_excel('WieringermeerData_Meteo.xlsx', index_col=0, parse_dates=['datetime'])








'''Filtering Data to Compute only one year'''


filtered = meteo.query('20180101 <= index <= 20191231')
meteo_P     = np.array(filtered['rain_station'])
meteo_Epot  = np.array(filtered['pEV'])
meteo_T     = np.array(filtered['temp'])

tfinal = len(meteo_P)

#Surface water flux conditions/forcings:

def BndWTop(tOut, meteo_P):  
    bndW = -1 * meteo_P[np.ceil(tOut).astype(int)-1]
    return bndW 

def BndETop(tOut, bPar):
    bndE = meteo_Epot[np.ceil(tOut).astype(int)-1] * bPar.cropfactor
    return bndE











'''Root Uptake Functions'''



def alpha(hw, sPar):
    h1 = sPar.h1
    h2 = sPar.h2
    h3 = sPar.h3
    h4 = sPar.h4
    
    
 
    ah = 0 * (np.logical_or(hw<=h4, hw>=h1)) \
            +1*(np.logical_and(hw>=h3, h2>=hw)) +\
                ((hw-h4)/(h3-h4))*(np.logical_and(h3>hw, h4<hw)) + \
                    ((hw-h1)/(h2-h1))*(np.logical_and(h1>hw, h2<hw))
    return ah

def Lrv(mdim, sPar):
    zN = mDim.zN
    zeta = sPar.zetaroot
    varrho = sPar.varrho
    L_rv = (zeta *np.exp(varrho * zN))*(zN<=0) + 0*(zN>0)

    return L_rv

def S(tOut, hw, sPar, mDim):
    nIN = mDim.nIN
    nN = mDim.nN
    S = np.ones(([nIN, hw.shape[1]]), dtype=hw.dtype)
    ii = np.arange(0, nN)
    root_length = Lrv(mDim, sPar)
    dz = mDim.dzN[1]
    a = alpha(hw, sPar)
    beta = np.zeros(np.shape(zN), dtype=hw.dtype)
    beta[ii] = root_length[ii] * dz / np.sum(root_length * dz)
    
    S = a * beta * (meteo_Epot[np.ceil(tOut).astype(int)-1]*bPar.cropfactor)
    
    return S

def qr(Theta_w, Theta_res, Lrv):
    


    num = 2.63e-3 * np.exp(62*(Theta_w - Theta_res))
    denom = 6.68 - np.log(Lrv)
    
    qr = num/denom
    qr = qr*(qr<=0) + qr*1e18*(qr>0)
            
    return qr










'''Unsaturated Flow'''


def Seff_calc(hw,sPar):
    hc = -hw
    Seff = (1 + ((hc > 0)*hc*sPar.a) ** sPar.n) ** (-sPar.m)
    return Seff 

def theta_w_calc(hw, sPar):
    theta_sat   = sPar.theta_sat
    theta_res   = sPar.theta_res
    Seff        = Seff_calc(hw,sPar)
    thetaw      = Seff*(theta_sat-theta_res)+theta_res
    return thetaw

def CCmplxF(hw, sPar):   #complex derivative, article

    dh           = np.sqrt(np.finfo(float).eps)
    hcmplx       = hw.real + 1j * dh
    theta        = theta_w_calc(hcmplx, sPar)
    C            = theta.imag / dh   
    return C

def Cprime_calc(hw,sPar):
    Sw  = Sw_calc(hw,sPar)
    C   = CCmplxF(hw,sPar)
    theta_sat = sPar.theta_sat
    beta_w  = 1.5e-10
    rho_w   = 1000
    g       = 9.81
    cv      = sPar.Cv
    Ssw     = rho_w * g * (cv + theta_sat * beta_w)
    Cprime  = C + Sw * Ssw
    return Cprime

def krw_calc(hw,sPar):
    Seff = Seff_calc(hw,sPar)
    krw=Seff**3
    return krw

def kIN_calc(hw,sPar,mDim):
   
    #import model dimensions
    nr,nc = hw.shape
    nIN = mDim.nIN
    
    #import parameters
    ksat = sPar.ksat
    krw = krw_calc(hw,sPar)
    kN = ksat*krw
    
    kIN = np.zeros([nIN,nc], dtype=hw.dtype)
    kIN[0]=kN[0]
    ii = np.arange(1, nIN - 1)    
    kIN[ii] = np.minimum(kN[ii - 1], kN[ii])
    kIN[nIN - 1] = kN[nIN - 2]
    return kIN

def Sw_calc(hw, sPar):
    thetaw = theta_w_calc(hw,sPar)
    theta_sat = sPar.theta_sat
    Sw = thetaw/theta_sat
    return Sw

def WaterFlux(t, hw, sPar, mDim, bPar):
    
    #import model dimensions
    #print(hw.shape)
    nr,nc=hw.shape
    nIN = mDim.nIN
    dzN = mDim.dzN
 
    
    #initialize arrays
    qw = np.zeros((nIN,nc),dtype=hw.dtype)

    ii = np.arange(1,nIN-1)
  
    #calculate variables:
    kIN = kIN_calc(hw, sPar, mDim)
    qw[ii] = -kIN[ii]*((hw[ii]-hw[ii-1])/dzN[ii-1] + 1)
    
    #water flux at top layer (given)
    qw[nIN-1] = BndWTop(t, meteo_P)

    
    #Bottom layer Robin condition
    #qw[0] = -kIN[0] #fist condition to test, it's like this because hw=0
    qw[0] = -bPar.ksatBot*(hw[0] - bPar.hwBot) #second contiditon to test
    
    return qw

def DivWaterFlux(t,hw,sPar,mDim,bPar):
    
    #importing model dimensions:
    nr,nc=hw.shape
    nN = mDim.nN
    
    #initializing arrays
    divhw = np.zeros([nN,nc], dtype=hw.dtype)
    
    #Calculating variables
    Cprime = Cprime_calc(hw,sPar)

    #Calculating water flux
    qw = WaterFlux(t,hw,sPar,mDim,bPar)
    ii = np.arange(0,nN)
    
    #Calculate divergence of flux
    divhw[ii] = -(qw[ii+1]-qw[ii])/(dzIN[ii] * Cprime[ii])

    return divhw

def Integrate(tRange, iniSt, sPar, mDim, bPar):
    
    def dYdt(t, hw):
        if len(hw.shape)==1:
            hw = hw.reshape(mDim.nN,1)
        
        if scenario[run_scenario] == scenario[0]:
            rates = DivWaterFlux(t, hw, sPar, mDim, bPar)
        
        else:
            rates = (DivWaterFlux(t, hw, sPar, mDim, bPar) - S(t, hw, sPar, mDim)) / Cprime_calc(hw, sPar) 
        
        return rates
        
        
    
    def jacFun(t,y):
        if len(y.shape)==1:
            y = y.reshape(mDim.nN,1)
    
        nr, nc = y.shape
        dh = np.sqrt(np.finfo(float).eps)
        ycmplx = y.copy().astype(complex)
        ycmplx = np.repeat(ycmplx,nr,axis=1)
        c_ex = np.ones([nr,1])* 1j*dh
        ycmplx = ycmplx + np.diagflat(c_ex,0)
        dfdy = dYdt(t, ycmplx).imag/dh
        return sp.coo_matrix(dfdy)
    
    # solve rate equatio
    t_span = [tRange[0],tRange[-1]]
    int_result = spi.solve_ivp(dYdt, t_span, iniSt.squeeze(), 
                                t_eval=tRange, 
                                method='BDF', vectorized=True, jac=jacFun, 
                                rtol=1e-4) 
    
    
    return int_result
    
    
















'''Main''' 

# Domain
nIN = 101
# soil profile until 1m meter depth
zIN = np.linspace(-5, 0, num=nIN).reshape(nIN, 1)  #initial soil profile from 15m depth to 0, one column and 151 rows
# nIN = np.shape(zIN)[0]
zN = np.zeros(nIN - 1).reshape(nIN - 1, 1) #empty soil profile with 1 column and 150 rows
zN[0, 0] = zIN[0, 0] #first value of soil profile is equal to first element of initial soil profile
zN[1:nIN - 2, 0] = (zIN[1:nIN - 2, 0] + zIN[2:nIN - 1, 0]) / 2 #creating staggered grid
zN[nIN - 2, 0] = zIN[nIN - 1] #final value is the same as initial profile
nN = np.shape(zN)[0] #lenght of soil profile (100)
ii = np.arange(0, nN - 1)
dzN = (zN[ii + 1, 0] - zN[ii, 0]).reshape(nN - 1, 1)
dzIN = (zIN[1:, 0] - zIN[0:-1, 0]).reshape(nIN - 1, 1)

# collect model dimensions in a pandas series: mDim
mDim = {'zN' : zN,
        'zIN' : zIN,
        'dzN' : dzN,
        'dzIN' : dzIN,
        'nN' : nN,
        'nIN' : nIN
        }
mDim = pd.Series(mDim)

    
# Parameters


if scenario[run_scenario] == scenario[0]:
    print(f'running: {scenario[run_scenario]}')
    
    #Van Gnuchten 1980s Params
    
    a = 2           #1/m
    n = 3      #
    m = 1-1/n       #[-]
    eps = 0.95
    
    sPar = {'a': a, 
            'n': n,
            'm': m,
            'ksat': 6.2, #m/d 
            'theta_sat': eps,
            'theta_res': 0.04,
            'eps': eps,
            'Cv': 2.6e-7 
            }

elif scenario[run_scenario] == scenario[1]:
    print(f'running: {scenario[run_scenario]}')
    a = 2           #1/m
    n = 3      #
    m = 1-1/n       #[-]
    eps = 0.95
    sPar = {'h1':-0.1,                      #water heads for root uptake 
            'h2':-0.8,
            'h3':-9,
            'h4':-15,
            'a': a, 
            'n': n,
            'm': m,
            'ksat': 6.2, #m/d 
            'theta_sat': eps,
            'eps': eps, 
            'theta_res': 0.04,
            'zetaroot': 7500,
            'varrho': 20.15,
            'Cv': 2.6e-7 
            }
elif scenario[run_scenario] == scenario[2]:
    print(f'running: {scenario[run_scenario]}')
    a = 2           #1/m
    n = 3      #
    m = 1-1/n       #[-]
    eps = 0.95
    sPar = {'h1':-0.1,                      #water heads for root uptake 
            'h2':-0.8,
            'h3':-9,
            'h4':-15,
            'a': a, 
            'n': n,
            'm': m,
            'ksat': 6.2, #m/d 
            'theta_sat': eps,
            'eps': eps,
            'theta_res': 0.04,
            'zetaroot': 7500,
            'varrho': 2,
            'Cv': 2.6e-7 
            }

elif scenario[run_scenario] == scenario[3]:
    print(f'running: {scenario[run_scenario]}')
    a = 2           #1/m
    n = 3      #
    m = 1-1/n       #[-]
    eps = 0.41
    sPar = {'h1':-0.1,                      #water heads for root uptake 
            'h2':-0.8,
            'h3':-9,
            'h4':-15,
            'a': a, 
            'n': n,
            'm': m,
            'ksat': 305, #m/d 
            'theta_sat': eps,
            'theta_res': 0.06,
            'zetaroot': 7500,
            'varrho': 2,
            'Cv': 1e-8 
            }

else:
    print('Error:  no valid scenario chosen in line 26, enter a value from 1 to 3.')
    sys.exit()

sPar = pd.Series(sPar)

bPar = {'hwBot': -1, #m external head
        'ksatBot':0.01, #1/d
        'cropfactor':0.9
        }  

bPar = pd.Series(bPar)

# Initial Conditions
wl = -1 #depth water table
hwIni = wl - zN

# Time Discretization
tOut = np.logspace(-14, np.log10(tfinal), num=tfinal)       
mt.tic()
int_result = Integrate(tOut, hwIni, sPar, mDim, bPar) 
mt.toc()

tOut = int_result.t  # time
nOut = np.shape(tOut)[0]


hw = int_result.y



# Dirichlet boundary condition: write boundary temperature to output.
if int_result.success:
    print('Integration has been successful')

qw = WaterFlux(tOut, hw, sPar, mDim, bPar)
thetaw = theta_w_calc(hw, sPar)
if scenario[run_scenario] != scenario[0]:
    qr = qr(thetaw, sPar.theta_res, Lrv(mDim, sPar))
    S = S(tOut, hw, sPar, mDim)


plt.close('all')
fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10, 6.18))
for i in range(0, nIN-1, 5):
    ax1.plot(tOut, hw[i,:])# label=f'depth -{zN[i][0]:.3f}m')
ax1.plot(tOut, hw[-1,:], '--', label=f'depth -{zN[-1][0]}m')
fig1.suptitle(f'Pressure Head: {scenario[run_scenario]}')
ax1.set_title('Water pressure head (ODE)')
ax1.set_xlabel('time (days)')
ax1.set_ylabel('Water pressure head [m]')
ax1.axvline(0, label='01/01/2018', ls='-.')
ax1.axvline(365, label='01/01/2019', ls='-.')
ax1.axvline((365*2)-1, label='31/12/2019', ls='-.')
ax1.legend(loc='best')
fig1.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/Pressure_Head_{scenario[run_scenario]}.png', dpi=350, format='png')

fig11, ax11 = plt.subplots(nrows=1, ncols=1, figsize=(13, 6.18))
ax11.plot(meteo_P, label='Precipitation in m')
ax11.plot(meteo_Epot, label='Potential evaporation in m')
ax11.legend(loc=0)
ax11.set_ylim(0, 0.05)
ax11.set_xlim(tOut[0], tOut[-1])
ax11.axvline(0, label='01/01/2018', ls='-.')
ax11.axvline(365, label='01/01/2019', ls='-.')
ax11.axvline((365*2)-1, label='31/12/2019', ls='-.')
ax11.legend(loc='best')
fig11.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/P_EPOT_{scenario[run_scenario]}.png', dpi=350, format='png')


fig2, ax2 = plt.subplots(figsize=(10, 10))
for ii in np.arange(0, nOut, 5):
    ax2.plot(hw[:, ii], zN, '-')
ax2.plot(hw[:, 0], zN, linestyle='dashed', label=f'day:{0}m')
fig2.suptitle(f'Pressure head: {scenario[run_scenario]}')
ax2.set_title(r'$h_w$ vs. depth (ODE)')
ax2.set_ylabel('depth [m]')
ax2.set_xlabel('Water pressure head [m]')
ax2.legend(loc='best')
fig2.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/Pressure head_{scenario[run_scenario]}.png', dpi=350, format='png')


'''Water FLux'''


fig3, ax3 = plt.subplots(figsize=(10, 10))
for ii in np.arange(3, nOut, 5):
    day = ii
    ax3.plot(qw[:, ii], zIN, '-')
ax3.plot(qw[:, 0], zIN, linestyle='dashed', label=f'day:{0}')
ax3.set_title('Water Flux vs. depth (ODE)')
ax3.set_ylabel('depth [m]')
ax3.set_xlabel('Water flux [m/d]')
ax3.legend(loc='best')
fig3.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/WaterFlux_{scenario[run_scenario]}.png', dpi=350, format='png')


'''Water Content'''



fig4, ax4 = plt.subplots(figsize=(10, 10))
for ii in np.arange(nOut-1, 10, -1):
    ax4.plot(thetaw[:, ii], zN, '-')
ax4.plot(thetaw[:, 0], zN, '-', label=f'day:{0}m')
fig4.suptitle(f'Volumetric Water Content: {scenario[run_scenario]}')
ax4.set_title(r'Water content ($\Theta_{w}$) vs. depth (ODE)')
ax4.set_ylabel('depth [m]')
ax4.set_xlabel(r'$\Theta_w [m/d]$')
ax4.legend(loc='best')
fig4.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/thetaw_{scenario[run_scenario]}.png', dpi=350, format='png')

if scenario[run_scenario] != scenario[0]:
    '''Sink'''
    
    fig5, ax5 = plt.subplots(1, 1, figsize=(10, 6.18), sharex=True)
    for i in range(1, nIN-1, 1):
        ax5.plot(tOut, -S[i,:])
    fig5.suptitle(f'Sink: {scenario[run_scenario]}')
    ax5.plot(tOut, S[1,:], '--', label=f'depth -{zN[-1][0]:.3f}m')
    ax5.set_title(f'Sink, varrho = {sPar.varrho}')
    ax5.set_xlabel('time (days)')
    ax5.set_ylabel(r'Sink [$m/d$] ')
    ax5.axvline(0, label='01/01/2018', ls='-.')
    ax5.axvline(365, label='01/01/2019', ls='-.')
    ax5.axvline((365*2)-1, label='31/12/2019', ls='-.')
    ax5.legend(loc='best')
    fig5.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/sink_{scenario[run_scenario]}.png', dpi=350, format='png')
    
        



'''     __  Root Water Plotting   __  '''

if scenario[run_scenario] != scenario[0]:

    '''Alpha z'''
    
    plt.figure(figsize=(10, 10))
    plt.plot(alpha(hw, sPar), zN)
    plt.title('Water Uptake Stress Function')
    plt.suptitle(f'Alpha: {scenario[run_scenario]}')
    plt.xlabel(r'$alpha[-]$')
    plt.ylabel(r'$z[m]$')
    plt.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/alphaz_{scenario[run_scenario]}.png', dpi=350, format='png')
    
    
    '''head alpha'''
    head_array = [sPar.h1, sPar.h2, sPar.h3, sPar.h4]
    head_y_array = [0, 1, 1, 0] 
    xticks = [-20 ,-15, -10, -5, 0]
    plt.figure(figsize=(10, 6.18))
    plt.plot(head_array, head_y_array)
    plt.scatter(head_array, head_y_array, marker='x' , color='red')
    plt.title('Water Uptake Stress Function')
    plt.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/headalpha{scenario[run_scenario]}.png', dpi=350, format='png')
    plt.xlabel(r'$h_w[m]$')
    plt.ylabel(r'$ \alpha [-]$')
    plt.xticks(xticks);

    '''Lrv'''
    plt.figure(figsize=(10, 10))
    plt.plot(Lrv(mDim, sPar), zN)
    plt.suptitle(f'Root Density: {scenario[run_scenario]}')
    plt.title(r'$L_{rv} / RD   [m/m3]$')
    plt.ylim(-5, 0)
    plt.xlabel(r'$L_{rv}[m/m^3]$')
    plt.ylabel(r'$ z[-]$');
    plt.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/LRV_{scenario[run_scenario]}.png', dpi=350, format='png')
    
    
    '''qroot'''
    
    
    plt.figure(figsize=(10, 10))
    for i in range(0, len(tOut)):
        plt.plot((-1e-5/qr[:, i]), zN, label=f'day: {i}')
    
    plt.suptitle(f'Water Uptake Rate: {scenario[run_scenario]}')
    plt.title(r'Root water flux $q_r$')
    plt.ylim(-1.2, 0)
    # plt.xlim(-0.01, 1)
    plt.xlabel(r'$q_{r}[m^3 m^{-1} d^{-1}]$')
    plt.ylabel(r'$ z[-]$')
    plt.savefig(f'C:/Users/creeb/Documents/TUDelft/Q4/CIE4365_Modelling_Coupled_Processes/Module 4 - Personal/Coding/outputs/{scenario[run_scenario]}/q_root_{scenario[run_scenario]}.png', dpi=350, format='png');











    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    