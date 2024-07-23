
import numpy as np
from numpy import sin, cos, pi, sqrt
import numdifftools as nd

class Pulse:
    def __init__(self,tf,omega0,tao,sigma, param5=None):

        
        self.tf = tf
        self.omega0 = omega0
        self.tao = tao
        self.sigma = sigma
        self.t = np.linspace(0, self.tf, 1000)
        self.param5 = param5


        
    def Omega1(self):
        return sqrt(self.Omega2()**2 + self.Omega3()**2) 
 
    def Omega2(self):

        return self.omega0 * np.exp(-((self.t-self.tao-self.tf/2)/self.sigma)**2)

    def Omega3(self):

        return self.omega0 * np.exp(-((self.t+self.tao-self.tf/2)/self.sigma)**2)

    def Omega4(self):
   
        return sqrt(self.Omega2()**2 + self.Omega3()**2) 

    def Omega2_t(self):

        def function(t):
            return self.omega0 * np.exp(-((t-self.tao-self.tf/2)/self.sigma)**2)
        return function
    def Omega3_t(self):
        
        def function(t):
            return self.omega0 * np.exp(-((t+self.tao-self.tf/2)/self.sigma)**2)
        return function
    
# 定义导数函数
    def dOmega2_dt(self):
        return nd.Derivative(self.Omega2_t(),step=1e-9)

    def dOmega3_dt(self,):
        return nd.Derivative(self.Omega3_t(),step=1e-9)

    def Omegacd(self):
        dOmega2_dt = self.dOmega2_dt()
        dOmega3_dt = self.dOmega3_dt()
        return (dOmega2_dt(self.t)*self.Omega3() - dOmega3_dt(self.t)*self.Omega2())/(self.Omega2()**2 + self.Omega3()**2)
        