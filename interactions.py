"""Functions for handling particle collisions and decays"""
import numpy as np
from constants import pi
from particle import Particle


def lorentzBoost(vector,beta=0,direction=[0,0,0]):
    """Calculate the lorentz boosted vector for given beta and direction of boost"""
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
    theta = np.random.random_sample()*2*pi #decay angle
    m0 = particle.mass
    m1 = products[0].mass
    m2 = products[1].mass
    pmag1 = np.sqrt((m0**2+m1**2-m2**2)**2 - 4*m0**2*m1**2) / (2*m0)
    p1 = [pmag1*np.cos(theta),0,pmag1*np.sin(theta)]
    p2 = [-pmag1*np.cos(theta),0,-pmag1*np.sin(theta)]

    # Boost momenta to lab frame
    fourmom1 = lorentzBoost([products[0].energy]+p1, particle.beta,particle.dir)
    fourmom2 = lorentzBoost([products[1].energy]+p2, particle.beta,particle.dir)
    products[0].momentum = fourmom1[1:]
    products[1].momentum = fourmom2[1:]

    return products
