"""Code to generate hadron shower from primary"""
from random import random,choice
from constants import pi
from particle import Particle


def randomValue(start,stop=None):
    """Returns a random value in the range [start,stop)"""
    if stop is None:
        stop = start
        start = 0
    return start+random()*(stop-start)

def generateRandomPrimary():
    """Returns a randomized primary particle"""
    # primaryTypes = ["proton"]
    particleType = "proton"
    position = [0,0,1000]
    energy = 1000
    theta = pi
    phi = 0

    return Particle(particleType,pos=position,E=energy,theta=theta,phi=phi)
