"""Code to generate hadron shower from primary"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from constants import pi
from MCmethods import randomDistance, chooseEnergy
from particle import Particle
from atmosphere import density, getAtmosphericNucleus, getAirCrossSection, getCollisionInverseCDF
from interactions import decay, collision, getDecayInverseCDF



def generatePrimary(**kwargs):
    """Returns a randomized primary particle with optional set attributes"""
    # Random values
    particleType = "proton"
    height = 500000
    position = [0,0,height]
    ke = chooseEnergy(minimum=1000)
    theta = pi
    phi = 0

    # Set values
    for key,val in kwargs.items():
        if key=="type":
            particleType = val
        elif key=="height":
            height = val
            position = [0,0,height]
        elif key=="minE":
            ke = chooseEnergy(minimum=val)
            if "energy" in kwargs.keys():
                ke = kwargs["energy"]
        elif key=="energy":
            ke = val
        elif key=="theta":
            theta = val
        elif key=="phi":
            phi = val
        else:
            print("Warning: argument",key,"in generatePrimary not recognized")

    return Particle(particleType,pos=position,KE=ke,theta=theta,phi=phi)


def getNextInteraction(particle):
    """Return a propagation length for the particle and the particle it
    interacts with after the propagation (decay returns "decay" as target,
    continued propagation returns None as target)"""
    if particle.lifetime is not None:
        decayInvCDF = getDecayInverseCDF(particle.lifetime,particle.beta)
        decayLength = randomDistance(decayInvCDF)
    else:
        decayLength = None
    if not("mu" in particle.type):
        sigma = getAirCrossSection(particle.energy)
        z = particle.position[2]
        theta = particle.theta
        collisionInvCDF = getCollisionInverseCDF(sigma,z,theta)
        collisionLength = randomDistance(collisionInvCDF)
    else:
        collisionLength = None

    if decayLength is None:
        return collisionLength, getAtmosphericNucleus(particle.position)
    elif collisionLength is None:
        return decayLength, "decay"
    else:
        if decayLength<collisionLength:
            return decayLength, "decay"
        else:
            return collisionLength, getAtmosphericNucleus(particle.position)


def propagate(particle,ceiling=None):
    """Propagate the particle and return particle that caused it to stop
    (decay returns "decay")"""
    distance, target = getNextInteraction(particle)
    if particle.position[2]+distance*particle.direction[2]<0:
        distance = -1*particle.position[2]/particle.direction[2]+.1
        target = None
    if ceiling is not None and \
       particle.position[2]+distance*particle.direction[2]>ceiling:
        distance = particle.position[2]/particle.direction[2]+.1
        target = None
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


def generateShower(primary,drawShower=False,maxIterations=1000):
    """Generates a full hadron shower and returns any muons that reach the surface"""
    #Setup
    particles = [primary]
    finished = False
    propagationParticles = ["pi+","pi-","mu+","mu-","p+","n0"]
    if drawShower:
        vertices = {primary.id: [[x for x in primary.position]]}
        colors = {primary.id: drawColor(primary.type)}
        markers = {primary.id: drawMarker(primary.type)}

    # Set a ceiling above which particles can be assumed to escape
    ceiling = 2*primary.position[2]

    # Loop until all propagating particles reach the ground
    loopCount = 0
    while not(finished):
        loopCount += 1
        products = []
        # Propagate and interact any particles that need it
        for particle in particles:
            if particle.type in propagationParticles and \
               (particle.position[2]>0 and particle.position[2]<ceiling):
                target = propagate(particle,ceiling)
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
        #     print(particle.id,particle.type,particle.position[2])

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
            # print("Stopped after",loopCount,"iterations")
            break

    # Get the muons
    muons = []
    for particle in particles:
        if (particle.type=="mu+" or particle.type=="mu-") and \
           particle.position[2]<0:
            muons.append(particle)

    # Plot the shower development
    if drawShower:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for particleId,points in vertices.items():
            xvals = [pos[0] for pos in points]
            yvals = [pos[1] for pos in points]
            zvals = [pos[2] for pos in points]
            ax.plot(xs=xvals,ys=yvals,zs=zvals,
                    color=colors[particleId],marker=markers[particleId])
        # Only show plot below first interaction point
        plotHeight = vertices[0][1][2]*1.2
        ax.set_zbound(-10,plotHeight)
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
    proton = generatePrimary()
    mus = generateShower(proton,drawShower=True)
    print(len(mus),"muons reached the ground")
