# -*- coding: utf-8 -*-
"""
current_models - library of ionic current models implemented in Python

Created on Mon Apr 10 16:30:04 2017

@author: Oliver Britton
"""

import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

" Abstract classes "
 
class Gate(object):
    def __init__(self):
        pass
    
class Current(object):
    def __init__(self, initial_conditions=None):
        self.initial_conditions = initial_conditions # Add initial conditions here
        pass

" Voltage clamp generator functions "


" //--Nav models--\\ "

" -- Nav 1.7 models -- "

def nav17vw(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Human Nav 1.7 from Vasylyev Waxman "
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    
    alpha_m = 10.22 - 10.22/(1 + np.exp((v+7.19)/15.43)) # Rate for closed -> open (sort of)
    beta_m = 23.76/(1 + np.exp((v+70.37)/14.53)) # Rate for open->closed
    
    """
    Parameters for kinetics - rate constant (3), 2 voltage shifts, 2 slope coefficients.
    """

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    alpha_h = 0.0744/(1 + np.exp((v+99.76)/11.07))
    beta_h = 2.54 - 2.54/(1 + np.exp((v+7.8)/10.68))

    hinf = alpha_h/(alpha_h + beta_h)
    htau = 1/(alpha_h + beta_h)
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    return [dm, dh]

def nav17cw(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Rat? Nav 1.7 from Choi Waxman 2011 "
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    s = Y[2]
    
    alpha_m = 15.5/(1 + np.exp(-(v-5)/(12.08)))
    beta_m = 35.2/(1 + np.exp((v+72.7)/16.7))

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    alpha_h = 0.38685/(1 + np.exp((v+122.35)/15.29))
    beta_h = -0.00283 + 2.00283/(1 + np.exp(-(v+5.5266)/12.70195)) # Rate is negative if v = -inf?

    hinf = alpha_h/(alpha_h + beta_h)
    htau = 1/(alpha_h + beta_h)

    alpha_s = 0.00003 + 0.00092/(1 + np.exp((v+93.9)/16.6))
    beta_s = 132.05 - 132.05/(1 + np.exp((v-384.9)/28.5))

    sinf = alpha_s/(alpha_s + beta_s)
    stau = 1/(alpha_s + beta_s)
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    ds = (sinf-s)/stau
    
    return [dm, dh, ds]
    
" -- Nav 1.8 models -- "
def nav18hw(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Human Nav 1.8 from Huang Waxman 20(14?) "
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    
    alpha_m = 7.35 - 7.35/(1 + np.exp((v+1.38)/10.9))
    beta_m = 5.97/(1 + np.exp((v+56.43)/18.26))

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    alpha_h = 0.011 + 1.39/(1 + np.exp((v+78.04)/11.32))
    beta_h = 0.56 - 0.56/(1 + np.exp((v-21.82)/20.03))

    hinf = alpha_h/(alpha_h + beta_h)
    htau = 1/(alpha_h + beta_h)

    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    
    return [dm, dh]

def nav18tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Rat? Nav 1.8 used in Tigerholm model "
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    s = Y[2]
    u = Y[3]
    
    alpha_m = 2.85 - 2.839/(1 + np.exp((v-1.159)/13.95))
    beta_m = 7.6205/(1 + np.exp((v+46.463)/8.8289))

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    hinf = 1/(1+np.exp((v+32.2)/4))
    htau = 1.218 + 42.043*np.exp(-((v+38.1)**2)/(2*15.19**2))

    alpha_s = 0.001 * 5.4203 / (1 + np.exp((v+79.816)/16.269))
    beta_s  = 0.001 * 5.0757 / (1 + np.exp(-(v+15.968)/11.542))

    sinf = 1/(1+np.exp((v+45.0)/8))
    stau = 1/(alpha_s + beta_s)

    alpha_u =  0.002 * 2.0434 / (1 + np.exp((v+67.499)/19.51))
    beta_u =  0.002 * 1.9952 / (1 + np.exp(-(v+30.963)/14.792))

    uinf = 1/(1+np.exp((v+51.0)/8))
    utau = 1.0/(alpha_u + beta_u) 
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    ds = (sinf-s)/stau
    du = (uinf-u)/utau
    
    return [dm, dh, ds, du]
    
def nav18cw(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Nav 1.8 model used in Choi Waxman 2011 "
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    
    alpha_m = 2.85 - 2.839/(1 + np.exp((v-1.159)/13.95))
    beta_m = 7.6205/(1 + np.exp((v+46.463)/8.8289))

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    hinf = 1/(1+np.exp((v+32.2)/4))
    htau = 1.218 + 42.043*np.exp(-((v+38.1)**2)/(2*15.19**2))
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    
    return [dm, dh]
    
" -- Nav 1.9 models -- "

def nav19hw(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Nav 1.9 model from Huang Waxman 2014"
    m = Y[0]
    h = Y[1]
    s = Y[2]
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    alpha_m = 0.751/(1 + np.exp(-(v+32.26)/13.71))
    beta_m = 5.68/(1 + np.exp((v+123.71)/13.94))
    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    alpha_h = 0.082/(1 + np.exp((v+113.69)/17.4))
    beta_h = 0.24/(1 + np.exp(-(v-10.1)/17.2))
    hinf = alpha_h/(alpha_h + beta_h)
    htau = 1/(alpha_h + beta_h)
    
    alpha_s = 0.019/(1 + np.exp((v+154.51)/11.46))
    beta_s = 0.000376/(1 + np.exp(-(v+60.92)/15.79))
    sinf = alpha_s/(alpha_s + beta_s)
    stau = 1/(alpha_s + beta_s)
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    ds = (sinf-s)/stau
        
    return [dm, dh, ds]
    
def nav19md(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Nav 1.9 model from Maingret 2008"
    m = Y[0]
    h = Y[1]
    s = Y[2]
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    alpgha
    
    return [dm, dh, ds]

" Kv models "

def kdr_tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Tigerholm version of the Sheets et al. IKdr model "
    " Model was developed from data recorded at 21 oC "
    
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    n = Y[0]
    q10 = 1.0#3.3 # Preserved in case it is useful but disabled
    
    if v > -31.0:
        tau = 0.16+0.8*np.exp(-0.0267*(v+11))
    else:
        tau = 1000*(0.000688 + 1/(np.exp((v+75.2)/6.5) + np.exp(-(v-131.5)/(34.8))))
		
    ninf = 1/(1 + np.exp(-(v+45)/15.4))
    ntau = tau/q10
    
    dn = (ninf-n)/ntau
    return [dn]
    
def km_tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    """ Tigerholm version of the IM current. Current is from multiple sources:
    The voltage dependence of steady-state activation forthe KM current is from
    Maingret et al. (2008), which was derived from Passmore 2003. The KM channel activation has a fast and a slow 
    time constant as described by Passmore et al. (2003). To account for the 
    two time constants, weimplemented one fast (nf) and one slow (ns) gate, 
    combined as follows.
    """
    # g = gbar * (0.25*ns + 0.75*nf)
    v = voltage_clamp_func(t,voltage_clamp_params)
    ns = Y[0]
    nf = Y[1]
    q10 = 1.0#3.3 # Preserved in case it is useful but disabled
    
    if v < -60.0:
        nstau = 219.0*q10
    else:
        nstau = 13.0*v + 1000.0*q10
        
    nftau_alpha = 0.00395*np.exp((v+30.0)/40.0)
    nftau_beta = 0.00395*np.exp(-(v+30.0)/20.0)*q10
    nftau = 1.0/(nftau_alpha + nftau_beta)
    
    ninf = 1.0/(1.0 + np.exp(-(v+30.0)/6.0)) # Threshold is around -30 mV
    
    dns = (ninf-ns)/nstau
    dnf = (ninf-nf)/nftau
    
    return [dns,dnf]
    
def ka_tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    """ Tigerholm version of IA.
    """
    # g = gbar * n * h
    v = voltage_clamp_func(t,voltage_clamp_params)
    n = Y[0]
    h = Y[1]
    q10 = 1.0#3.3 # Preserved in case it is useful but disabled
    
    ninf = (1.0/(1.0 + np.exp(-(v+5.4+15)/16.4)))**4
    ntau = 0.25 + 10.04*np.exp((-(v+24.67)**2)/(2*34.8**2))*q10
		
    hinf = 1.0/(1.0 + np.exp((v+49.9 + 15.0)/4.6))
    htau = 20.0 + 50.0 * np.exp((-(v+40.0)**2)/(2.0*40.0**2))*q10
    
    # Trap for htau following Sheets /ChoiWaxman/Tigerholm - set it to 5 ms if less than 5 ms
    if htau < 5.0:
        htau = 5.0

    dn = (ninf-n)/ntau
    dh = (hinf-h)/htau
    
    return [dn,dh]

" HCN models "
def hcn_kn(Y,t,voltage_clamp_func,voltage_clamp_params):
    """ 
    Kouranova Ih model with non-specific current (reversal potential should be set at -30 mV 
    """

    v = voltage_clamp_func(t,voltage_clamp_params)
    n_s = Y[0]
    n_f = Y[1]

    ninf_s = 1/(1 + np.exp((v+87.2)/9.7))
    ninf_f = ninf_s

    if v > -70.0:
        tau_ns = 300.0 + 542.0 * np.exp((v+25.0)/20.0)
        tau_nf = 140.0 + 50.0 * np.exp(-(v+25.0)/20.0)
    else:
        tau_ns = 2500.0 + 100.0 * np.exp((v+240.0)/50.0)
        tau_nf = 250.0 + 12.0 * np.exp((v+240.0)/50.0)

    dns = (ninf_s - n_s)/tau_ns
    dnf = (ninf_f - n_f)/tau_nf

    return [dns, dnf]
    
def hcn_tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    """
    Tigerholm version of the Kouranova Ih model which is identical except
    that when you calculate the current you don't use a nonspecific reversal potential and instead split the current between Na+ and K+, 50/50.    
    """
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    n_s = Y[0]
    n_f = Y[1]
    
    ninf_s = 1/(1 + np.exp((v+87.2)/9.7))
    ninf_f = ninf_s

    if v > -70.0:
        tau_ns = 300.0 + 542.0 * np.exp((v+25.0)/20.0)
        tau_nf = 140.0 + 50.0 * np.exp(-(v+25.0)/20.0)
    else:
        tau_ns = 2500.0 + 100.0 * np.exp((v+240.0)/50.0)
        tau_nf = 250.0 + 12.0 * np.exp((v+240.0)/50.0)

    dns = (ninf_s - n_s)/tau_ns
    dnf = (ninf_f - n_f)/tau_nf
    
    return [dns, dnf]

"""
    # ena, ek, + or -?
    Ih_na = 0.5 * g_h (0.5*n_s + 0.5*n_f) * (Vm + ena)
    Ih_k = 0.5 * g_h * (0.5*n_s + 0.5*n_f) * (Vm + ek) 

"""

" Test models "
def nav17test(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Human Nav 1.7 from Vasylyev Waxman "
    
    v = voltage_clamp_func(t,voltage_clamp_params)
    
    m = Y[0]
    h = Y[1]
    
    alpha_m = 10.22 - 10.22/(1 + np.exp((v+7.19)/15.43)) # Rate for closed -> open (sort of)
    beta_m = 23.76/(1 + np.exp((v+70.37)/14.53)) # Rate for open->closed
    
    """
    Parameters for kinetics - rate constant (3), 2 voltage shifts, 2 slope coefficients.
    """

    minf = alpha_m/(alpha_m + beta_m)
    mtau = 1/(alpha_m + beta_m)

    alpha_h = 0.0744/(1 + np.exp((v+99.76)/11.07))
    beta_h = 2.54 - 2.54/(1 + np.exp((v+7.8)/10.68))

    hinf = alpha_h/(alpha_h + beta_h)
    htau = 1/(alpha_h + beta_h)
    
    dm = (minf-m)/mtau
    dh = (hinf-h)/htau
    return [dm, dh]
    

 