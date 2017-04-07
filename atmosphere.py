"""Atmospheric model code"""
from particle import Particle
from random import random

def density(height):
    """Returns the number density (in m^-3) of the atmosphere at a given height"""
    return 0.02504e27

def getAtmosphericNucleus():
    """Returns an atomic nucleus based on the atmospheric composition"""
    atomSeed = random()
    if atomSeed>1-.78:
        return Particle("nitrogen")
    elif atomSeed>.01:
        return Particle("oxygen")
    else:
        return Particle("argon")
