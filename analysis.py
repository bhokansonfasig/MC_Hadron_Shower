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

    titleString = "Muon Lateral Distribution\nfor "+info["count"]
    if info["isotropic"]:
        titleString += " Isotropic"
    titleString += " Events with "
    if info["energyType"]=="minE":
        titleString += "E>"
    else:
        titleString += "E="
    titleString += info["energyValue"]

    logmin = np.log10(np.min(muonDistances))
    logmax = np.log10(np.max(muonDistances))

    plt.hist(muonDistances,log=True,bins=np.logspace(logmin,logmax))
    plt.gca().set_xscale("log")
    plt.title(titleString)
    plt.xlabel("Distance to shower core (m)")
    plt.ylabel("Flux (arbitrary units)")
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


def plotEnergyDistribution(dataFileName,plotName=None):
    """Generates histogram of energies of muons reaching the ground"""
    pFile = open(dataFileName,'rb')
    showerResults = pickle.load(pFile)
    info = interpretFilename(dataFileName)

    muonEnergies = []
    for i,muons in enumerate(showerResults):
        for muon in muons:
            if muon.ke>50:
                muonEnergies.append(muon.energy)

    titleString = "Muon Energy Distribution\nfor "+info["count"]
    if info["isotropic"]:
        titleString += " Isotropic"
    titleString += " Events with "
    if info["energyType"]=="minE":
        titleString += "E>"
    else:
        titleString += "E="
    titleString += info["energyValue"]

    logmin = np.log10(np.min(muonEnergies))
    logmax = np.log10(np.max(muonEnergies))

    plt.hist(muonEnergies,log=True,bins=np.logspace(logmin,logmax))
    plt.gca().set_xscale("log")
    plt.title(titleString)
    plt.xlabel("Muon energy (MeV)")
    plt.ylabel("Flux (arbitrary units)")
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


def plotMomentumDistribution(dataFileName,plotName=None):
    """Generates histogram of momenta of muons reaching the ground"""
    pFile = open(dataFileName,'rb')
    showerResults = pickle.load(pFile)
    info = interpretFilename(dataFileName)

    muonMomenta = []
    for i,muons in enumerate(showerResults):
        for muon in muons:
            if muon.ke>50:
                muonMomenta.append(muon.Pmag)

    titleString = "Muon Momentum Distribution\nfor "+info["count"]
    if info["isotropic"]:
        titleString += " Isotropic"
    titleString += " Events with "
    if info["energyType"]=="minE":
        titleString += "E>"
    else:
        titleString += "E="
    titleString += info["energyValue"]

    logmin = np.log10(np.min(muonMomenta))
    logmax = np.log10(np.max(muonMomenta))

    plt.hist(muonMomenta,log=True,bins=np.logspace(logmin,logmax))
    plt.gca().set_xscale("log")
    plt.title(titleString)
    plt.xlabel("Muon momentum (MeV)")
    plt.ylabel("Flux (arbitrary units)")
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


def plotPrimaryEnergies(num,minE=100,plotName=None):
    """Generates histogram of energies of proton primaries"""
    scaledE, letter = scaleValue(minE,1e6)
    energyString = "E>"+str(int(scaledE))+letter+"eV"

    primaryEnergies = np.zeros(num)
    for i in range(num):
        proton = generatePrimary(minE=minE)
        primaryEnergies[i] = proton.energy

    titleString = "Primary Energy Distribution\nfor "+str(num)
    titleString += " Events with "
    titleString += energyString

    logmin = np.log10(primaryEnergies.min())
    logmax = np.log10(primaryEnergies.max())

    plt.hist(primaryEnergies,log=True,bins=np.logspace(logmin,logmax))
    plt.gca().set_xscale("log")
    plt.title(titleString)
    plt.xlabel("Proton energy (MeV)")
    plt.ylabel("Flux (arbitrary units)")
    if plotName is not None:
        if not(isinstance(plotName,str)):
            plotName = titleString.replace(" ","_")
        plt.savefig(plotName)
    plt.show()


if __name__ == '__main__':
    # plotSingleShower(energy=1e9)

    # count = 10000
    # for energy in [100,1e3,1e6,1e9,1e12]:
    #     generateDataset(count,minE=energy)

    # plotPrimaryEnergies(100000,plotName="100000_minE100MeV_primaries.png")

    fileBase = "10000_setE1PeV"
    plotNumberCounts(fileBase+".pickle",plotName=fileBase+"_numcts.png")
    plotLateralDistribution(fileBase+".pickle",plotName=fileBase+"_latdist.png")
    plotEnergyDistribution(fileBase+".pickle",plotName=fileBase+"_energies.png")
    plotMomentumDistribution(fileBase+".pickle",plotName=fileBase+"_momenta.png")
