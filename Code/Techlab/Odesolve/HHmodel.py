# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 16:46:52 2016

@author: Oliver Britton
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def hodgkinhuxley(y, t, params):

    V = y[0]
    m = y[1]
    h = y[2]
    n = y[3]
    
    Vm = -60
    C = 1
    
    gNa = 120  
    eNa = 115
    
    alpha_m = (2.5-0.1*(V-Vm))/(exp(2.5-0.1*(V-Vm)) -1)
    beta_m = 4*exp(-(V-Vm)/18)
    dm = alpha_m * (1-m) - beta_m * m
    
    alpha_h = 0.07*exp(-(V-Vm)/20)
    beta_h = 1./(exp(3.0-0.1*(V-Vm))+1)
    dh = alpha_h * (1-h) - beta_h * h
     
    INa = gNa * m^3 * h * ((V-Vm)-eNa)

    gK = 36
    eK = -12

    alpha_n = (0.1-0.01*(V-Vm)) ./ (exp(1-0.1*(V-Vm)) -1)
    beta_n = 0.125*exp(-(V-Vm)/80)
    dn = alpha_n * (1-n) - beta_n * n
    
    IK =  gK * n^4 * ((V-Vm)-eK)
    
    gL=0.3 
    eL = 10.6
    
    Ileak = gL * ((V-Vm)-eL)
    
    Istim = 0
    
    derivs = 