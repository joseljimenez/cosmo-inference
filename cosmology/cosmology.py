# Here all the cosmology
from numbers import Integral
import numpy as np
from scipy import integrate
#import logging

C_LIGHT = 299792.458 #km/s

class Cosmology:
    def __init__(self, Ob0, H0, As, alpha):
        self.Ob0 = Ob0
        self.H0 = H0
        self.As = As
        self.alpha = alpha
    
    def f_unified(self, z):
        return (self.As + (1 - self.As)*(z + 1)**(3*(1 + self.alpha)))**(1/(1 + self.alpha))

    def H(self, z):
        return (self.H0)*np.sqrt(self.Ob0*(1+z)**3 + (1-self.Ob0)*self.f_unified(z))

    def integrand(self, z):
        return 1/self.H(z)

    def DA(self, z):
        integral, _ = integrate.quad(self.integrand, 0, z)
        return (C_LIGHT/(1 + z))*integral
    
    def DM(self, z):
        return (1 + z)*self.DA(z)

    def DH(self, z):
        return C_LIGHT/(self.H(z))

    def Dv(self, z):
        return (z * self.DM(z)**2 * self.DH(z))**(1 / 3)

    def Dl(self, z):
        return (1 + z)**2 * self.DA(z)
    
    def mu(self, z):
        return 5 * np.log10(self.Dl(z)) + 25

    def mb(self, z, M):
        return self.mu(z) + M
