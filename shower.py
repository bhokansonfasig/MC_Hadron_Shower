"""Code to generate hadron shower from primary"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from constants import pi
from particle import Particle
from atmosphere import density, getAtmosphericNucleus
from interactions import decay, collision


def randomValue(start,stop=None):
    """Returns a random value in the range [start,stop)"""
    if stop is None:
        stop = start
        start = 0
    return start+np.random.random_sample()*(stop-start)

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
    interacts with after the propagation (decay returns "decay" as target,
    continued propagation returns None as target)"""
    interactionLength = 100 #m
    if "mu" in particle.type or "pi" in particle.type:
        target = "decay"
    else:
        if particle.ke<140*5:
            target = None
        else:
            target = getAtmosphericNucleus(particle.position)
    return interactionLength, target

def propagate(particle):
    """Propagate the particle and return particle that caused it to stop
    (decay returns "decay")"""
    distance, target = getNextInteraction(particle)
    for i in range(len(particle.position)):
        particle.position[i] += distance * particle.direction[i]
    return target


def interact(particle,target=None):
    """Interact the particle with the target (or "decay") and return the
    products. If target is None, do nothing and return the particle"""
    if target is None:
        return [particle]
    elif target is "decay":
        decay(particle)
    elif target.type=="N" or target.type=="O" or target.type=="Ar":
        collision(particle,target)
    else:
        print("Warning: no supported interaction between",particle.type,
              "and",target.type)
        return [particle,target]


def generateShower(drawShower=False):
    """Generates a full hadron shower and returns any muons that reach the surface"""
    #Setup
    primary = generateRandomPrimary()
    particles = [primary]
    finished = False
    propagationParticles = ["pi+","pi-","mu+","mu-","p+","n0"]
    if drawShower:
        vertices = {primary.id: [[x for x in primary.position]]}
        colors = {primary.id: drawColor(primary.type)}

    # Loop until all propagating particles reach the ground
    while not(finished):
        products = []
        # Propagate and interact any particles that need it
        for particle in particles:
            if particle.type in propagationParticles:
                target = propagate(particle)
                if drawShower:
                    vertices[particle.id].append([x for x in particle.position])
                products.extend(interact(particle,target))
            else:
                products.append(particle)
        particles = products
        finished = True
        # Print particles at each step
        print("---")
        for particle in particles:
            print(particle.id,particle.type)
        # Add positions to the drawing dictionary
        if drawShower:
            for particle in particles:
                try:
                    vertices[particle.id].append([x for x in particle.position])
                except KeyError:
                    vertices[particle.id] = [[x for x in particle.position]]
                    colors[particle.id] = drawColor(particle.type)
        # Check to see if all propagating particles have reached the ground
        for particle in particles:
            if particle.type in propagationParticles and particle.position[2]>0:
                finished = False
                break
    # Get the muons
    muons = []
    for particle in particles:
        if particle.type=="mu+" or particle.type=="mu-":
            muons.append(particle)
    # Plot the shower development
    if drawShower:
        # for pos in vertices[0]:
        #     print(pos)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for particleId,points in vertices.items():
            xvals = [pos[0] for pos in points]
            yvals = [pos[1] for pos in points]
            zvals = [pos[2] for pos in points]
            ax.plot(xs=xvals,ys=yvals,zs=zvals,color=colors[particleId])#,marker='.')
        plt.show()
    return muons


def drawColor(particleType):
    """Return the color the particle should be drawn in"""
    if particleType[:2]=="mu":
        return 'b'
    elif particleType[:2]=="pi":
        return 'g'
    elif particleType=="p+" or particleType=="n0":
        return 'r'
    elif particleType[:2]=="nu":
        return 'c'
    else:
        return 'k'



if __name__ == '__main__':
    generateShower(drawShower=True)
