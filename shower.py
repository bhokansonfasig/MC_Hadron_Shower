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


def getPropagationLength(particle):
    """Return a propagation length for the particle"""
    interactionLength = 100 #m
    return interactionLength

def propagate(particle):
    """Propagate the particle"""
    distance = getPropagationLength(particle)
    for i in range(len(particle.position)):
        particle.position[i] += distance * particle.direction[i]


def interact(particle,target=None):
    """Interact the particle with the target (or decay by default)
    and return the products"""
    if target is None:
        return [particle]
    else:
        return [particle,target]


def generateShower():
    """Generates a full hadron shower and returns any muons that reach the surface"""
    primary = generateRandomPrimary()
    particles = [primary]
    finished = False
    while not(finished):
        products = []
        for particle in particles:
            if particle.type=="p+" or particle.type=="n0":
                propagate(particle)
                products.extend(interact(particle,Particle("nitrogen")))
            elif particle.type=="pi+" or particle.type=="pi-":
                propagate(particle)
                products.extend(interact(particle,Particle("nitrogen")))
            elif particle.type=="mu+" or particle.type=="mu-":
                propagate(particle)
                products.extend(interact(particle))
            else:
                products.append(particle)
        particles = products
        finished = True
        print("---")
        for particle in particles:
            print(particle.type,particle.position)
        for particle in particles:
            propagationParticles = ["pi+","pi-","mu+","mu-","p+","n0"]
            if particle.type in propagationParticles and particle.position[2]>0:
                finished = False
                break


