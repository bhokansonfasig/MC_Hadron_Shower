"""Functions used to generate random numbers for Monte Carlo"""
from numpy.random import random_sample, normal
from numpy import sqrt, log, pi, sin, cos, arccos, dot


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


def isotropicAngles(downgoing=False):
    """Generates theta and phi with isotropic distribution"""
    if downgoing:
        costheta = -1*random_sample()
    else:
        costheta = random_sample()*2-1
    theta = arccos(costheta)
    phi = random_sample()*2*pi
    return theta,phi


def pointInCircle(radius=1):
    """Returns random points x,y within a circle of set radius"""
    r2 = random_sample()*radius**2
    theta = random_sample()*2*pi
    return sqrt(r2)*cos(theta),sqrt(r2)*sin(theta)


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


def chooseMultiplicity(labE,totalKE):
    """Choose a pion multiplicity value based on energies"""
    # Expected value of the multiplicity
    expected = 6*log(labE)-24
    # Maximum value of the multiplicity = KE/(3*pion mass)
    maximum = totalKE / (139.57018+139.57018+134.9766)

    # Choose from a gaussian distribution until multiplicity is a reasonable value
    mult = -1
    while mult<0 or mult>maximum:
        mult = int(normal(loc=expected,scale=sqrt(expected)))

    return mult


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
    seed = random_sample()
    return inverseCDF(seed)


def chooseEnergy(minimum=100):
    """Returns a random energy value based on an E^-2.7 power law spectrum
    with a minimum energy of 1 GeV"""
    seed = random_sample()
    normalization = minimum**1.7
    return ((1-seed)/normalization)**(-1/1.7)
