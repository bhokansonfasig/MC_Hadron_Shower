"""Atmospheric model code"""
import numpy as np
from MCmethods import random_sample
from particle import Particle

def density(height,scaleHeight=8000):
    """Returns the number density (in m^-3) of the atmosphere at a given height"""
    return 2.504e25*np.exp(-height/scaleHeight)

def getAtmosphericNucleus(position):
    """Returns an atomic nucleus based on the atmospheric composition"""
    atomSeed = random_sample()
    if atomSeed>=1-.78: #78% nitrogen
        return Particle("nitrogen",pos=position)
    elif atomSeed>=.01: #21% oxygen
        return Particle("oxygen",pos=position)
    else:              #1% argon
        return Particle("argon",pos=position)


def getCollisionInverseCDF(crossSection,particleHeight,particleTheta,scaleHeight=8000):
    """Returns a function that will give the distance value for a provided value
    of the cumulative density function"""
    def invCDF(val):
        """Calculates the inverse cumulative density function for the
        atmospheric model"""
        a = density(0)*crossSection
        b = scaleHeight/np.cos(particleTheta)
        c = np.exp(-particleHeight/scaleHeight)
        if b>0:
            n = 1-np.exp(-a*b*c)
        else:
            n = 1
        # Catch negative logarithms and assume they should be nearly log(0)
        if 1+np.log(1-n*val)/(a*b*c)<0:
            return 1e99

        return -b*np.log(1+ np.log(1-n*val)/(a*b*c))

    return invCDF


def getAirCrossSection(energy):
    """Returns the cross section of a particle with the air at given energy"""
    s = 2*energy/1000*13 #GeV^2

    sigmapp = 32.4 - 1.2*np.log(s) + 0.21*np.log(s)**2

    sigmapair = sigmapp*4+100

    #Convert from mb to m^2
    return sigmapair*1e-31
