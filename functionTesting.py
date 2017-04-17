"""Tests for functions across project"""
import numpy as np
import matplotlib.pyplot as plt
from particle import Particle
from interactions import lorentzBoost,decay,collision
from MCmethods import randomDistance
from atmosphere import getInverseCDF

def testDecay():
    """Test decay function with rest-frame charged pion decay"""
    print("Charged pion decay test-")
    pion = Particle("pi+")
    products = decay(pion)
    for prod in products:
        if prod.type=="mu+":
            print("  Muon kinetic energy:",prod.ke)
    print("       Expected value: 4.12")

def testLorentzBoost():
    """Test lorentzBoost function with a basic boost"""
    print("Lorentz boost test-")
    p = [5,0,0,2]
    p2 = lorentzBoost(p, .6,[0,0,-1])
    expected = [7.75,0,0,6.25]
    for i,x in enumerate(["E","x","y","z"]):
        print("  Boosted "+x+"-component:",p2[i])
        print("       Expected value:",expected[i])


def testCollision():
    """Test basic kinematics of proton-Nitrogen collision"""
    print("Basic collision test-")
    proton = Particle("p+",momentum=[0,0,-10000])
    nitrogen = Particle("N")
    e0 = proton.energy + nitrogen.energy
    p0 = proton.momentum
    products = collision(proton,nitrogen)
    print("  Products:",end=' ')
    for prod in products:
        print(prod.type,end=" ")
    print()
    print("  Initial energy:",e0)
    ef = 0
    pf = [0,0,0]
    for prod in products:
        ef += prod.energy
        for i in range(3):
            pf[i] += prod.momentum[i]
            pf[i] = round(pf[i],3)
    print("    Final energy:",ef)

    print("  Initial momentum:",p0)
    print("    Final momentum:",pf)


def testRandomDistance():
    """Test random distance distribution for atmospheric model"""
    trials = 10000
    nbins = 50

    sigma = 2e-25
    invLambda = 0.02504e27*sigma
    z = 200000
    theta = 2/3*np.pi

    distances = np.zeros(trials)
    function = getInverseCDF(sigma,z,theta)
    for i in range(trials):
        distances[i] = randomDistance(function)
    plt.hist(distances,bins=nbins,label="Data",normed=True)

    xMax = int(-z/np.cos(theta))
    xVals = np.arange(0,xMax,xMax/1000)
    expVals = np.exp(-(z+xVals*np.cos(theta))/8000)
    yVals = invLambda*expVals*np.exp(-invLambda*8000/np.cos(theta) * \
                                     (np.exp(-z/8000)-expVals))
    plt.plot(xVals,yVals,label="Theory")

    plt.legend()
    plt.show()


if __name__ == '__main__':
    # testDecay()
    # testLorentzBoost()
    # testCollision()
    testRandomDistance()
