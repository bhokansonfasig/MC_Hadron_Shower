"""Functions for doing analysis of simulated hadron showers"""
import pickle
import datetime
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


def interpretFilename(filename):
    """Returns values of primary based on filename"""
    extensionIndex = filename.rfind(".p")
    filenameBase = filename[:extensionIndex]
    bits = filenameBase.split("_")

    info = {"count": bits[0],
            "energyType": bits[1][:4],
            "energyValue": bits[1][4:],
            "theta": None, "phi": None,
            "isotropic": False}
    if len(bits)==3:
        info["isotropic"] = True
    elif len(bits)==4:
        info["theta"] = bits[2][5:]
        info["phi"] = bits[3][3:]
    elif len(bits)==5:
        info["theta"] = bits[2][5:]
        info["phi"] = bits[3][3:]
        info["isotropic"] = True

    return info


def plotSingleShower(plotName=None,**kwargs):
    """Plots one shower in 3D, kwargs passed on to generatePrimary"""
    proton = generatePrimary(**kwargs)
    muons = generateShower(proton,drawShower=True,plotName=plotName)
    print(len(muons),"muons reached the ground")


def generateDataset(num,minE=100,setE=None,isotropic=False,theta=None,phi=None):
    """Generate num showers and return muons from each shower"""
    if setE is None:
        scaledE, letter = scaleValue(minE,1e6)
        energyString = "minE"+str(int(scaledE))+letter+"eV"
    else:
        scaledE, letter = scaleValue(setE,1e6)
        energyString = "setE"+str(int(scaledE))+letter+"eV"
    print(datetime.datetime.now().strftime("%H:%M"),"- Generating",end=" ")
    if isotropic:
        print("isotropic",end=" ")
    else:
        print("set",end=" ")
    print("dataset with",energyString[:4]+"="+energyString[4:])

    showerResults = []
    tracker = 10
    for i in range(num):
        if setE is None:
            proton = generatePrimary(minE=minE,isotropic=isotropic,
                                     theta=theta,phi=phi)
        else:
            proton = generatePrimary(energy=setE,isotropic=isotropic,
                                     theta=theta,phi=phi)
        muons = generateShower(proton)
        if 100*i/num>=tracker:
            print("      -",str(tracker)+"%","@",
                  datetime.datetime.now().strftime("%H:%M"))
            tracker += 10
        showerResults.append(muons)

    filename = str(num)+"_"
    filename += energyString
    if theta is not None and phi is not None:
        filename += "_theta"+str(round(theta,2))+"_phi"+str(round(phi,2))
    if isotropic:
        filename += "_isotropic"
    filename += ".pickle"

    print("      - Saving to",filename)

    pFile = open(filename,'wb')
    pickle.dump(showerResults,pFile,-1)
    pFile.close()

    return filename


def plotNumberCounts(dataFileName,plotName=None):
    """Generates histogram of number of muons reaching the ground"""
    pFile = open(dataFileName,'rb')
    showerResults = pickle.load(pFile)
    info = interpretFilename(dataFileName)

    muonCounts = np.zeros(int(info["count"]))
    for i,muons in enumerate(showerResults):
        for muon in muons:
            if muon.ke>50:
                muonCounts[i] += 1

    print(np.count_nonzero(muonCounts),"events with muons")

    titleString = "Muon Number Distribution\nfor "+info["count"]
    if info["isotropic"]:
        titleString += " Isotropic"
    titleString += " Events with "
    if info["energyType"]=="minE":
        titleString += "E>"
    else:
        titleString += "E="
    titleString += info["energyValue"]

    plt.hist(muonCounts,np.arange(-.5,muonCounts.max()+1.5,1))
    plt.title(titleString)
    plt.xlabel("Number of muons reaching the ground")
    plt.xlim([-.5,muonCounts.max()+.5])
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


def plotLateralDistribution(dataFileName,rmax=None,plotName=None):
    """Generates histogram of radii of muons reaching the ground"""
    pFile = open(dataFileName,'rb')
    showerResults = pickle.load(pFile)
    info = interpretFilename(dataFileName)

    muonDistances = []
    for i,muons in enumerate(showerResults):
        for muon in muons:
            if muon.ke>50:
                r = np.sqrt(muon.position[0]**2+muon.position[1]**2)
                if rmax is None or r<rmax:
                    muonDistances.append(r)

    titleString = "Muon Number Distribution\nfor "+info["count"]
    if info["isotropic"]:
        titleString += " Isotropic"
    titleString += " Events with "
    if info["energyType"]=="minE":
        titleString += "E>"
    else:
        titleString += "E="
    titleString += info["energyValue"]

    plt.hist(muonDistances,bins=50)
    plt.title(titleString)
    plt.xlabel("Distance to shower core (m)")
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


if __name__ == '__main__':
    # plotSingleShower(energy=1e9)
    # muonNumberCounts(1000,minE=200)
    # lateralDistribution(1000,setE=1e12,rmax=2000)

    count = 10000
    for energy in [1e3,1e6,1e9,1e12]:
        generateDataset(count,setE=energy)
