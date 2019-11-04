import numpy as np
import scipy.integrate

from Methods.currents import get_voltage_clamp_vector, voltage_clamp

"""
Classes that implement current models with all parameters named so we can manipulate them easily 
e.g. for sensitivity analyses

TO DOs:

* Possibly make a virtual base class if it is useful, but current models can be quite different from one another 
so I don't know if this will save time. If we repeat ourselves 3 times then I will do this.

* Add Q10 modelling, by coding in the model's base temperature, and creating Q10, temperature, and flag parameters to enable Q10 scaling. 

"""

class current(object):
        return
    

class nav17vw(object):
    " Nav 1.7 current model from Vasylev Waxman "
    
    def __init__(self, g = 0.18, e_rev=88.5):
    
        self.names = {'gates':['m', 'h'], 'parameters':['rates', 'vshifts', 'vwidths']}
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
        self.num_state_vars = len(self.initial_conditions)
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
            
    def calc_IV(self, v):
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
        
        
class kdr_tf(object):
    " IKdr model used in Tigerholm et al., 2103 from Sheets et al. 2007"
        
    def __init__(self, g = 0.0018, e_rev=-100.05):
        
        self.names = {'gates':['n'], 'parameters':['vshifts', 'vwidths']}
        self.vshifts = {'n':45.0}
        self.vwidths
        
        self.initial_conditions = [0]
        self.num_state_vars = len(self.initial_conditions)
        self.g = g
        self.e_rev = e_rev