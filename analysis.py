"""Functions for doing analysis of simulated hadron showers"""
import numpy as np
import matplotlib.pyplot as plt
from shower import generatePrimary, generateShower


def scaleValue(value,base=1):
    """Returns scaled value and corresponding letter"""
    if value>=1:
        factor = 1/1000
        letters = ["","k","M","G","T","P","E","Z","Y"]
    else:
        factor = 1000
        letters = ["m","mu","n","p","f","a","z","y"]
        value *= 1000
    value *= base
    valString = str(float(value))
    if 'e' not in valString and value>1:
        decimalIndex = valString.find('.')
        index = int((decimalIndex-1)/3)
        return value*factor**index, letters[index]
    else:
        index = 0
        while value<1 or value>=1000:
            index += 1
            value *= factor
        return value, letters[index]



def plotSingleShower(plotName=None,**kwargs):
    """Plots one shower in 3D, kwargs passed on to generatePrimary"""
    proton = generatePrimary(**kwargs)
    muons = generateShower(proton,drawShower=True,plotName=plotName)
    print(len(muons),"muons reached the ground")


def muonNumberCounts(num,minE=100,setE=None,plotName=None):
    """Generates histogram of number of muons reaching the ground"""
    muonCounts = np.zeros(num)
    for i in range(num):
        if setE is None:
            proton = generatePrimary(minE=minE)
        else:
            proton = generatePrimary(energy=setE)
        muons = generateShower(proton)
        for muon in muons:
            if muon.ke>50:
                muonCounts[i] += 1

    if setE is None:
        scaledE, letter = scaleValue(minE,1e6)
        energyString = "E>"+str(int(scaledE))+" "+letter+"eV"
    else:
        scaledE, letter = scaleValue(setE,1e6)
        energyString = "E="+str(int(scaledE))+" "+letter+"eV"

    plt.hist(muonCounts,np.arange(-.5,muonCounts.max()+1.5,1))
    plt.title("Muon Number Distribution\nfor "+str(num)+" Events"+\
              " with "+energyString)
    plt.xlabel("Number of muons reaching the ground")
    plt.xlim([-.5,muonCounts.max()+.5])
    if plotName is not None:
        plt.savefig(plotName)
    plt.show()


def lateralDistribution(num,minE=100,setE=None,rmax=None,plotName=None):
    """Generates histogram of radii of muons reaching the ground"""
    muonDistances = []
    for _ in range(num):
        if setE is None:
            proton = generatePrimary(minE=minE)
        else:
            proton = generatePrimary(energy=setE)
        muons = generateShower(proton)
        for muon in muons:
            if muon.ke>50:
                r = np.sqrt(muon.position[0]**2+muon.position[1]**2)
                if rmax is None or r<rmax:
                    muonDistances.append(r)

    if setE is None:
        scaledE, letter = scaleValue(minE,1e6)
        energyString = "E>"+str(int(scaledE))+" "+letter+"eV"
    else:
        scaledE, letter = scaleValue(setE,1e6)
        energyString = "E="+str(int(scaledE))+" "+letter+"eV"

    plt.hist(muonDistances,bins=50)
    plt.title("Muon Lateral Distribution\nfor "+str(num)+" Events"+\
              " with "+energyString)
    plt.xlabel("Distance to shower core (m)")
    if plotName is not None:
        plt.savefig(plotName)
    plt.show()


if __name__ == '__main__':
    # plotSingleShower(energy=1e9)
    # muonNumberCounts(100,setE=1000000)
    lateralDistribution(100,setE=1e6,rmax=None)
