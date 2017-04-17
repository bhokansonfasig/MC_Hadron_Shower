"""Functions for handling particle collisions and decays"""
import numpy as np
from constants import pi
from particle import Particle


def lorentzBoost(vector,beta=0,direction=[0,0,0]):
    """Calculate the lorentz boosted 4-vector for given beta and direction of boost"""
    nx = direction[0]
    ny = direction[1]
    nz = direction[2]
    gamma = 1/np.sqrt(1-beta**2)
    boost = [[gamma, -gamma*beta*nx, -gamma*beta*ny, -gamma*beta*nz],
             [-gamma*beta*nx, 1+(gamma-1)*nx**2, (gamma-1)*nx*ny, (gamma-1)*nx*nz],
             [-gamma*beta*ny, (gamma-1)*ny*nx, 1+(gamma-1)*ny**2, (gamma-1)*ny*nz],
             [-gamma*beta*nz, (gamma-1)*nz*nx, (gamma-1)*nz*ny, 1+(gamma-1)*nz**2]]
    return np.dot(boost,vector)


def rotate3D(vector,theta=0,phi=0):
    """Rotate 3-vector by spherical angles theta and phi"""
    rx = [[1,0,0],
          [0,np.cos(theta),-np.sin(theta)],
          [0,np.sin(theta),np.cos(theta)]]
    rz = [[np.cos(phi),-np.sin(phi),0],
          [np.sin(phi),np.cos(phi),0],
          [0,0,1]]
    return np.dot(rz,np.dot(rx,vector))


def randomWithSum(n,total):
    """Generates n uniformly distributed random numbers whose sum is total"""
    seeds = np.random.random_sample(n-1)
    seeds = [0,1]+list(seeds)
    seeds.sort()
    vals = []
    for i in range(n):
        vals.append(total*(seeds[i+1]-seeds[i]))
    return vals


def decay(particle):
    """Calculate kinematics of two body decay and return resulting particles.
    If the decay should have three particles, any neutrinos are ignored and
    the resulting one- or two-body decay is performed"""
    # Determine products
    if particle.type=="pi+":
        products = [Particle("mu+",pos=particle.position),
                    Particle("nuMu",pos=particle.position)]
    elif particle.type=="pi-":
        products = [Particle("mu-",pos=particle.position),
                    Particle("nuMuBar",pos=particle.position)]
    elif particle.type=="F-16":
        prodcuts = [Particle("oxygen-15",pos=particle.position),
                    Particle("proton",pos=particle.position)]
    # The following product lists ignore neutrinos
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
    elif particle.type=="C-14" or particle.type=="O-14":
        products = [Particle("nitrogen-14",pos=particle.position),
                    Particle("e-",pos=particle.position)]
    elif particle.type=="N-16":
        products = [Particle("oxygen-16",pos=particle.position),
                    Particle("e-",pos=particle.position)]
    elif particle.type=="O-15":
        products = [Particle("nitrogen-15",pos=particle.position),
                    Particle("e+",pos=particle.position)]
    elif particle.type=="Cl-40":
        products = [Particle("argon-40",pos=particle.position),
                    Particle("e-",pos=particle.position)]
    elif particle.type=="K-40":
        products = [Particle("calcium-40",pos=particle.position),
                    Particle("e-",pos=particle.position)]
    else:
        print("Warning:",particle.type,"decay not known")
        return [particle]

    # Kinematics in rest frame
    phi = np.random.random_sample()*2*pi #decay angle
    costheta = np.random.random_sample()*2-1 #rotation around interaction axis
    theta = np.arccos(costheta)
    m0 = particle.mass
    m1 = products[0].mass
    m2 = products[1].mass
    pmag1 = np.sqrt((m0**2+m1**2-m2**2)**2 - 4*m0**2*m1**2) / (2*m0)
    p1 = [pmag1*np.sin(theta)*np.cos(phi),
          pmag1*np.sin(theta)*np.sin(phi),
          pmag1*np.cos(theta)]
    p2 = [-p1[0],-p1[1],-p1[2]]

    # Boost momenta to lab frame
    products[0].momentum = p1
    products[1].momentum = p2
    beta = particle.beta
    direction = [-1*x for x in particle.dir]
    fourmom1 = lorentzBoost([products[0].energy]+p1, beta,direction)
    fourmom2 = lorentzBoost([products[1].energy]+p2, beta,direction)
    products[0].momentum = fourmom1[1:]
    products[1].momentum = fourmom2[1:]

    return products



def collision(particle,target):
    """Calculate kinematics of collision between particle and target.
    Returns products of the collision"""
    # Determine product types
    pionSeed = np.random.random_sample()
    if pionSeed>=2/3:   #1/3 pi+
        pionType = "pi+"
        if target.type=="N":
            productType = "carbon-14"
        elif target.type=="O":
            productType = "nitrogen-16"
        else:
            productType = "chlorine-40"
    elif pionSeed>=1/3: #1/3 pi-
        pionType = "pi-"
        if target.type=="N":
            productType = "oxygen-14"
        elif target.type=="O":
            productType = "fluorine-16"
        else:
            productType = "potassium-40"
    else:               #1/3 pi0
        pionType = "pi0"
        productType = target.type

    products = [Particle(particle.type,id=particle.id,pos=particle.position),
                Particle(productType,pos=particle.position),
                Particle(pionType,pos=particle.position)]

    # Determine kinematics by dividing non-mass energy randomly (in COM frame)
    # then determine the momentum vectors by forming a triangle with fixed side
    # lengths (assume first vector points along z-axis)

    # Boost to the COM frame before getting energies
    beta = particle.Pmag/(particle.energy+target.mass)
    direction = particle.dir
    fourmom1 = lorentzBoost([particle.energy]+particle.momentum, beta,direction)
    fourmom2 = lorentzBoost([target.energy]+target.momentum, beta,direction)

    productMassTotal = sum([particle.mass for particle in products])
    totalKE = fourmom1[0] + fourmom2[0] - productMassTotal
    if totalKE<=0:
        # Not enough energy for collision to do anything
        return [particle,target]
    # Loop through generation of momenta until they could possibly sum to zero
    realistic = False
    while not(realistic):
        kes = randomWithSum(3,totalKE)
        pmags = [0,0,0]
        for i,prod in enumerate(products):
            pmags[i] = np.sqrt(kes[i]**2 + 2*kes[i]*prod.mass)
        # Momenta could sum to zero if no one is greater than the sum of the others
        realistic = (pmags[0]<pmags[1]+pmags[2]) and \
                    (pmags[1]<pmags[0]+pmags[2]) and \
                    (pmags[2]<pmags[0]+pmags[1])

    # Set angles from x-axis (unrotated)
    xis = [0,
           pi-np.arccos((pmags[0]**2+pmags[1]**2-pmags[2]**2)/2/pmags[0]/pmags[1]),
           pi+np.arccos((pmags[0]**2+pmags[2]**2-pmags[1]**2)/2/pmags[0]/pmags[2])]

    # Rotate triangle of momenta by some theta and phi
    phi = np.random.random_sample()*2*pi #decay angle
    costheta = np.random.random_sample()*2-1 #rotation around interaction axis
    theta = np.arccos(costheta)
    sign = int(np.random.random_sample()<0.5)*2-1

    for i,mag in enumerate(pmags):
        unrotated = [mag*np.sin(xis[i]),0,mag*np.cos(xis[i])]
        products[i].momentum = list(sign*rotate3D(unrotated,theta,phi))

    # for i,x in enumerate("xyzE"):
    #     tot = 0
    #     for prod in products:
    #         if x=="E":
    #             tot += prod.energy
    #         else:
    #             tot += prod.momentum[i]
    #     print("  "+x+"-component total:",tot)

    # Boost momenta to the lab frame
    beta = particle.Pmag/(particle.energy+target.mass)
    direction = [-1*x for x in particle.dir]
    for i,prod in enumerate(products):
        fourmomentum = lorentzBoost([prod.energy]+prod.momentum, beta,direction)
        prod.momentum = fourmomentum[1:]

    return products

