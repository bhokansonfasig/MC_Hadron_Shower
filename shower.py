"""Code to generate hadron shower from primary"""
from random import random,choice
from constants import pi
from particle import Particle
from atmosphere import density, getAtmosphericNucleus


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
    energy = 10000
    theta = pi
    phi = 0

    return Particle(particleType,pos=position,E=energy,theta=theta,phi=phi)


def getNextInteraction(particle):
    """Return a propagation length for the particle and the particle it
    interacts with after the propagation (decay returns None as target)"""
    interactionLength = 100 #m
    if "mu" in particle.type or "pi" in particle.type:
        target = None
    else:
        if particle.ke<140*3/2:
            target = None
        else:
            target = getAtmosphericNucleus(particle.position)
    return interactionLength, target

def propagate(particle):
    """Propagate the particle and return particle that caused it to stop
    (decay returns None)"""
    distance, target = getNextInteraction(particle)
    for i in range(len(particle.position)):
        particle.position[i] += distance * particle.direction[i]
    return target


def interact(particle,target=None):
    """Interact the particle with the target (or decay by default)
    and return the products"""
    if target is None:
        if particle.type=="pi+":
            return [Particle("mu+",pos=particle.position,energy=particle.energy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("nuMu",pos=particle.position)]
        elif particle.type=="pi-":
            return [Particle("mu-",pos=particle.position,energy=particle.energy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("nuMuBar",pos=particle.position)]
        elif particle.type=="mu+":
            return [Particle("e+",pos=particle.position,energy=particle.energy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("nuE",pos=particle.position),
                    Particle("nuMuBar",pos=particle.position)]
        elif particle.type=="mu-":
            return [Particle("e-",pos=particle.position,energy=particle.energy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("nuEBar",pos=particle.position),
                    Particle("nuMu",pos=particle.position)]
        else:
            print("Warning: particle",particle.type,"does not decay")
            return [particle]
    elif target.type=="N" or target.type=="O" or target.type=="Ar":
        pionType = random()
        if pionType>=2/3:   #1/3 pi+
            splitEnergy = particle.energy/3
            return [Particle(particle.type,pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("n0",pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("pi+",pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi)]
        elif pionType>=1/3: #1/3 pi-
            splitEnergy = particle.energy/3
            return [Particle(particle.type,pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("p+",pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi),
                    Particle("pi-",pos=particle.position,energy=splitEnergy,
                             theta=particle.theta,phi=particle.phi)]
        else:               #1/3 pi0
            splitEnergy = particle.energy/3
            secondaryType = random()
            if secondaryType>=1/2: #1/2 proton
                return [Particle(particle.type,pos=particle.position,
                                 energy=splitEnergy,theta=particle.theta,
                                 phi=particle.phi),
                        Particle("p+",pos=particle.position,energy=splitEnergy,
                                 theta=particle.theta,phi=particle.phi),
                        Particle("pi0",pos=particle.position,energy=splitEnergy,
                                 theta=particle.theta,phi=particle.phi)]
            else:                  #1/2 neutron
                return [Particle(particle.type,pos=particle.position,
                                 energy=splitEnergy,theta=particle.theta,
                                 phi=particle.phi),
                        Particle("n0",pos=particle.position,energy=splitEnergy,
                                 theta=particle.theta,phi=particle.phi),
                        Particle("pi0",pos=particle.position,energy=splitEnergy,
                                 theta=particle.theta,phi=particle.phi)]
    else:
        print("Warning: no interaction between",particle.type,"and",target.type)
        return [particle,target]


def generateShower():
    """Generates a full hadron shower and returns any muons that reach the surface"""
    primary = generateRandomPrimary()
    particles = [primary]
    finished = False
    while not(finished):
        products = []
        for particle in particles:
            if particle.type=="p+" or particle.type=="n0" or \
               particle.type=="pi+" or particle.type=="pi-" or \
               particle.type=="mu+" or particle.type=="mu-":
                target = propagate(particle)
                products.extend(interact(particle,target))
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


