"""Atmospheric model code"""
import numpy as np
from MCmethods import random_sample
from particle import Particle

def density(height,scaleHeight=8000):
    """Returns the number density (in m^-3) of the atmosphere at a given height"""
    return 0.02504e27*np.exp(-height/scaleHeight)

def getAtmosphericNucleus(position):
    """Returns an atomic nucleus based on the atmospheric composition"""
    atomSeed = random_sample()
    if atomSeed>=1-.78: #78% nitrogen
        return Particle("nitrogen",pos=position)
    elif atomSeed>=.01: #21% oxygen
        return Particle("oxygen",pos=position)
    else:              #1% argon
        return Particle("argon",pos=position)


def getInverseCDF(crossSection,particleHeight,particleTheta,scaleHeight=8000):
    """Returns a function that will give the distance value for a provided value
    of the cumulative density function"""
    def invCDF(val):
        """Calculates the inverse cumulative density function for the
        atmospheric model"""
        a = density(0)*crossSection
        b = scaleHeight/np.cos(particleTheta)
        c = np.exp(-particleHeight/scaleHeight)
        return -b*np.log(1+ np.log(1-val)/(a*b*c))
    return invCDF
