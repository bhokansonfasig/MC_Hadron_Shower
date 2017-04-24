"""Functions for handling particle collisions and decays"""
import numpy as np
from constants import pi, c
from MCmethods import randomInRange, randomWithSum, isotropicAngles, randomMomentumTriangle, chooseMultiplicity
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
    theta,phi = isotropicAngles()
    # phi: decay angle, theta: rotation about interaction axis
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
    fourmom1 = lorentzBoost([products[0].energy]+list(p1), beta,direction)
    fourmom2 = lorentzBoost([products[1].energy]+list(p2), beta,direction)
    products[0].momentum = fourmom1[1:]
    products[1].momentum = fourmom2[1:]

    return products



def collision(particle,target):
    """Calculate kinematics of collision between particle and target.
    Returns products of the collision"""
    # Determine product types
    pionSeed = randomInRange(3)
    if pionSeed<1:   #1/3 pi+
        pionType = "pi+"
        if target.type=="N":
            productType = "carbon-14"
        elif target.type=="O":
            productType = "nitrogen-16"
        else:
            productType = "chlorine-40"
    elif pionSeed<2: #1/3 pi-
        pionType = "pi-"
        if target.type=="N":
            productType = "oxygen-14"
        elif target.type=="O":
            productType = "fluorine-16"
        else:
            productType = "potassium-40"
    else:            #1/3 pi0
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
    fourmom1 = lorentzBoost([particle.energy]+list(particle.momentum), beta,direction)
    fourmom2 = lorentzBoost([target.energy]+list(target.momentum), beta,direction)

    productMassTotal = sum([particle.mass for particle in products])
    totalKE = fourmom1[0] + fourmom2[0] - productMassTotal
    if totalKE<=0:
        # Not enough energy for collision to do anything
        return [particle,target]

    # Determine multiplicity and add that number of sets of pi+,pi-,pi0 to the
    # products
    mult = chooseMultiplicity(particle.energy,totalKE)
    setKEs = randomWithSum(mult+1,totalKE)

    # Get momenta for first set of products
    momenta = randomMomentumTriangle(setKEs[0],[prod.mass for prod in products])

    # Get momenta for all additional pions from multiplicity
    for i in range(mult):
        additionalPions = [Particle("pi+",pos=particle.position),
                           Particle("pi-",pos=particle.position),
                           Particle("pi0",pos=particle.position)]
        additionalMomenta = randomMomentumTriangle(setKEs[i+1],
                            [prod.mass for prod in additionalPions])
        products.extend(additionalPions)
        momenta.extend(additionalMomenta)

    # Assign the momenta to the particles
    for i,prod in enumerate(products):
        prod.momentum = momenta[i]

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
        fourmomentum = lorentzBoost([prod.energy]+list(prod.momentum), beta,direction)
        prod.momentum = fourmomentum[1:]

    return products



def getDecayInverseCDF(lifetime,beta):
    """Returns a function that will give the distance value for a provided value
    of the cumulative density function"""
    def invCDF(val):
        """Calculates the inverse cumulative density function for decay"""
        if beta**2==1:
            d = 1e300
        else:
            d = beta*c*lifetime/np.sqrt(1-beta**2)
        return -d*np.log(1-val)
    return invCDF
