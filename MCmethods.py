"""Functions used to generate random numbers for Monte Carlo"""
from numpy.random import random_sample, normal
from numpy import sqrt, pi, sin, cos, arccos, dot


def randomInRange(start,stop=None):
    """Returns a random value in the range [start,stop)"""
    if stop is None:
        stop = start
        start = 0
    return start+random_sample()*(stop-start)


def randomWithSum(n,total):
    """Generates n uniformly distributed random numbers whose sum is total"""
    seeds = random_sample(n-1)
    seeds = [0,1]+list(seeds)
    seeds.sort()
    vals = []
    for i in range(n):
        vals.append(total*(seeds[i+1]-seeds[i]))
    return vals


def isotropicAngles():
    """Generates theta and phi with isotropic distribution"""
    costheta = random_sample()*2-1
    theta = arccos(costheta)
    phi = random_sample()*2*pi
    return theta,phi


def randomMomentumTriangle(totalKE,masses):
    """Generates three momentum vectors in a random fashion such that
    their vector sum is zero"""
    assert(totalKE>0)
    assert(isinstance(masses,list))
    assert(len(masses)==3)
    # Loop through generation of momenta until they could possibly sum to zero
    realistic = False
    while not(realistic):
        kes = randomWithSum(3,totalKE)
        pmags = [0,0,0]
        for i,mass in enumerate(masses):
            pmags[i] = sqrt(kes[i]**2 + 2*kes[i]*mass)
        # Momenta could sum to zero if no one is greater than the sum of the others
        realistic = (pmags[0]<pmags[1]+pmags[2]) and \
                    (pmags[1]<pmags[0]+pmags[2]) and \
                    (pmags[2]<pmags[0]+pmags[1])

    # Set angles from x-axis (unrotated)
    xis = [0,
           pi-arccos((pmags[0]**2+pmags[1]**2-pmags[2]**2)/2/pmags[0]/pmags[1]),
           pi+arccos((pmags[0]**2+pmags[2]**2-pmags[1]**2)/2/pmags[0]/pmags[2])]

    # Rotate triangle of momenta by some theta and phi
    theta,phi = isotropicAngles()
    sign = int(randomInRange(2)<1)*2-1

    momenta = []
    for i,mag in enumerate(pmags):
        unrotated = [mag*sin(xis[i]),0,mag*cos(xis[i])]
        momenta.append(sign*rotate3D(unrotated,theta,phi))

    return momenta


def rotate3D(vector,theta=0,phi=0):
    """Rotate 3-vector by spherical angles theta and phi"""
    rx = [[1,0,0],
          [0,cos(theta),-sin(theta)],
          [0,sin(theta),cos(theta)]]
    rz = [[cos(phi),-sin(phi),0],
          [sin(phi),cos(phi),0],
          [0,0,1]]
    return dot(rz, dot(rx,vector))


def randomDistance(inverseCDF):
    """Returns a random distance for a particle to travel given some inverted
    cumulative density function of the distance (based on the atmospheric model)"""
    return inverseCDF(random_sample())
