import numpy as np
import scipy.integrate

from Methods.currents import get_voltage_clamp_vector, voltage_clamp

class nav17vw(object):
    
    def __init__(self, g = 0.18, e_rev=88.5):
    
        # Standard gate has two rates
        # Closed -> open (alpha)
        # Open -> closed (beta)
        self.rates = {  # Cannot be more than 1 rate as terms have to cancel
                                    'm':{'alpha':10.22, 'beta':23.76}, 
                                    'h':{'alpha':0.0744, 'beta':2.54} 
                                   }
        self.vshifts = {  
                                    'm':{'alpha':7.19, 'beta':70.37 }, 
                                    'h':{'alpha':99.76, 'beta':7.8 } 
                                   }
        self.vwidths  = {  
                                    'm':{'alpha':15.43, 'beta':14.53 }, 
                                    'h':{'alpha':11.07, 'beta':10.68 } 
                                   }
        # Initital conditions
        self.initial_conditions = [0,1]
        self.g = g
        self.e_rev = e_rev
    
    def calc_terms(self, v):
        rates = self.rates
        vshifts = self.vshifts
        vwidths = self.vwidths
        
        alpha_m = rates['m']['alpha'] - rates['m']['alpha']/(1 + np.exp((v + vshifts['m']['alpha'])/vwidths['m']['alpha']))
        beta_m = rates['m']['beta']/(1 + np.exp((v + vshifts['m']['beta'])/vwidths['m']['beta']))
        
        alpha_h = rates['h']['alpha']/(1 + np.exp((v + vshifts['h']['alpha'])/vwidths['h']['alpha']))
        beta_h = rates['h']['beta'] - rates['h']['beta']/(1 + np.exp((v + vshifts['h']['beta'])/vwidths['h']['beta']))
        
        minf = alpha_m/(alpha_m + beta_m)
        mtau = 1/(alpha_m + beta_m)

        hinf = alpha_h/(alpha_h + beta_h)
        htau = 1/(alpha_h + beta_h)

        return [minf, mtau, hinf, htau]
            
    def calc_I(self, v):
        " Calculate steady state current at a particular voltage "
        [minf, mtau, hinf, htau] = self.calc_terms(v) # Calc steady states
        I = self.g*minf*minf*minf*hinf*(v-self.e_rev)
        return I
        
    def calc_derivatives(self, Y, t, voltage_clamp_func, voltage_clamp_params):
        " ODE function to calculate state variables with an ODE solver"
        v = voltage_clamp_func(t,voltage_clamp_params)
        [minf, mtau, hinf, htau] = self.calc_terms(v)      
        m = Y[0]
        h = Y[1]
        dm = (minf-m)/mtau
        dh = (hinf-h)/htau
        return [dm, dh]
        
    def simulate(self, t, voltage_clamp_func, voltage_clamp_params):
        " Solve model to simulate current "
        v = get_voltage_clamp_vector(t,voltage_clamp_func, voltage_clamp_params)
        solution = scipy.integrate.odeint(self.calc_derivatives, self.initial_conditions, t, args=(voltage_clamp_func, voltage_clamp_params))
        m = solution[:,0]
        h = solution[:,1]
        I = self.g*m*m*m*h*(v-self.e_rev)
        return [v,I]