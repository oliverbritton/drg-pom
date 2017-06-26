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


" Kv models "

" HCN models "
def hcn_tf(Y,t,voltage_clamp_func,voltage_clamp_params):
    " Tigerholm version of the Kouranova Ih model (original Kouranova model is different - no separate Na and K parts) "
    
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
    

 