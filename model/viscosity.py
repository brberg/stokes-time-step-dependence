from __future__ import division
from dolfin import  *

def strain_rate(velocity):
    "Calculates the strain rate tensor from a velocity field."
    return 0.5*(nabla_grad(velocity) + nabla_grad(velocity).T)

class ViscosityCalculator:
    def __init__(self, creep_parameter):
        #Exponentiates given creep parameter based on Glen's flow law exponent n=3 for use in later calculation.
        self.creep_parameter_exponentiated = creep_parameter**(-1/3)

    def compute_viscosity(self, velocity):
        "Computes the viscosity based on glen's flow law and a given velocity field."
        viscosity_coefficient = 0.5*self.creep_parameter_exponentiated
        strainrate = strain_rate(velocity)
        viscosity = viscosity_coefficient/((0.5*strainrate[0,0]**2 + strainrate[0,1]**2 + 0.5*strainrate[1,1]**2 + 1E-10)**(1/3))
        return viscosity