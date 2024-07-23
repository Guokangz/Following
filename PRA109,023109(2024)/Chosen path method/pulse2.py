
import numpy as np
from numpy import sin, cos, pi, sqrt 
import numdifftools as nd

class Pulse2:
    def __init__(self,tf,gamma,delta,param5=None):

        self.tf = tf
        self.gamma = gamma
        self.delta = delta
        self.t = np.linspace(0,self.tf,10000)
        self.param5 = param5

    def Mu_t(self,t):
        return self.gamma/2*(1-cos(2*pi*t/self.tf))

    def theta_t(self,t):
        return pi/2+pi*t/(2*self.tf) -1/3*sin(2*pi*t/self.tf) +1/24*sin(4*pi*t/self.tf)

    
    def theta_dt(self):
        return nd.Derivative(lambda t :self.theta_t(t),step=1e-9)
    
    def Mu_dt(self):
        return nd.Derivative(lambda t :self.Mu_t(t),step=1e-9)
    
    def Omega_EP(self):
        theta_dt = self.theta_dt()(self.t)
        Mu_dt = self.Mu_dt()(self.t)
        theta_t = self.theta_t(self.t)
        Mu_t = self.Mu_t(self.t)
        return -theta_dt*sin(theta_t)*cos(Mu_t)/sin(Mu_t)-Mu_dt*cos(theta_t)
    
    def Omega_ES(self):
        theta_dt = self.theta_dt()(self.t)
        Mu_dt = self.Mu_dt()(self.t)
        theta_t = self.theta_t(self.t)
        Mu_t = self.Mu_t(self.t)
        return theta_dt*cos(theta_t)*cos(Mu_t)/sin(Mu_t)-Mu_dt*sin(theta_t)
    
    def Omega1(self):
        value = (self.delta**2*(self.Omega_EP()**2+self.Omega_ES()**2))**0.25
        value = np.nan_to_num(value, nan=0.0, posinf=0, neginf=0.0)
        return value
 
    def Omega2(self):

        value = abs(self.Omega_EP()*(self.delta**2/(self.Omega_EP()**2+self.Omega_ES()**2))**0.25)
        value = np.nan_to_num(value, nan=0.0, posinf=0, neginf=0.0)
        return value

    def Omega3(self):

        value =  abs(self.Omega_ES()*(self.delta**2/(self.Omega_EP()**2+self.Omega_ES()**2))**0.25)
        value = np.nan_to_num(value, nan=0.0, posinf=0, neginf=0.0)
        return value

    def Omega4(self):
        value = (self.delta**2*(self.Omega_EP()**2+self.Omega_ES()**2))**0.25
        value = np.nan_to_num(value, nan=0.0, posinf=0, neginf=0.0)
        return value
    

        