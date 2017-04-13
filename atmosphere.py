"""Atmospheric model code"""
import numpy as np
from particle import Particle

def density(height):
    """Returns the number density (in m^-3) of the atmosphere at a given height"""
    return 0.02504e27

def getAtmosphericNucleus(position):
    """Returns an atomic nucleus based on the atmospheric composition"""
    atomSeed = np.random.random_sample()
    if atomSeed>=1-.78: #78% nitrogen
        return Particle("nitrogen",pos=position)
    elif atomSeed>=.01: #21% oxygen
        return Particle("oxygen",pos=position)
    else:              #1% argon
        return Particle("argon",pos=position)
