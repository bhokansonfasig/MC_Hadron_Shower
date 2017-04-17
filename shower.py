"""Code to generate hadron shower from primary"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from constants import pi
from particle import Particle
from atmosphere import density, getAtmosphericNucleus
from interactions import decay, collision



def generateRandomPrimary():
    """Returns a randomized primary particle"""
    # primaryTypes = ["proton"]
    particleType = "proton"
    position = [0,0,1000]
    energy = 100000
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
        return decay(particle)
    elif target.type=="N" or target.type=="O" or target.type=="Ar":
        return collision(particle,target)
    else:
        print("Warning: no supported interaction between",particle.type,
              "and",target.type)
        return [particle,target]


def generateShower(drawShower=False,maxIterations=1000):
    """Generates a full hadron shower and returns any muons that reach the surface"""
    #Setup
    primary = generateRandomPrimary()
    particles = [primary]
    finished = False
    propagationParticles = ["pi+","pi-","mu+","mu-","p+","n0"]
    if drawShower:
        vertices = {primary.id: [[x for x in primary.position]]}
        colors = {primary.id: drawColor(primary.type)}
        markers = {primary.id: drawMarker(primary.type)}

    # Set a ceiling above which particles can be assumed to escape
    ceiling = primary.position[2]

    # Loop until all propagating particles reach the ground
    loopCount = 0
    while not(finished):
        loopCount += 1
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

        # # Print particles at each step
        # print("---")
        # for particle in particles:
        #     print(particle.id,particle.type)

        # Add positions to the drawing dictionary
        if drawShower:
            for particle in particles:
                try:
                    vertices[particle.id].append([x for x in particle.position])
                except KeyError:
                    vertices[particle.id] = [[x for x in particle.position]]
                    colors[particle.id] = drawColor(particle.type)
                    markers[particle.id] = drawMarker(particle.type)

        # Check to see if all propagating particles have reached the ground
        for particle in particles:
            if particle.type in propagationParticles and \
               (particle.position[2]>0 and particle.position[2]<ceiling):
                finished = False
                break
        if loopCount==maxIterations:
            print("Stopped after",loopCount,"iterations")
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
            ax.plot(xs=xvals,ys=yvals,zs=zvals,
                    color=colors[particleId],marker=markers[particleId])
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

def drawMarker(particleType):
    """Return the marker style the particle should be drawn with"""
    if particleType[:2]=="mu":
        return ','
    elif particleType[:2]=="pi":
        return ','
    elif particleType=="p+" or particleType=="n0":
        return ','
    elif particleType[:2]=="nu":
        return 'x'
    elif particleType.lower()!=particleType:
        return ','
    else:
        return 'x'



if __name__ == '__main__':
    generateShower(drawShower=True)
