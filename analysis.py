"""Functions for doing analysis of simulated hadron showers"""
import numpy as np
import matplotlib.pyplot as plt
from shower import generateRandomPrimary, generateShower


def plotSingleShower(primaryHeight):
    proton = generateRandomPrimary(500000)
    muons = generateShower(proton,drawShower=True)
    print(len(muons),"muons reached the ground")


def manyShowers(num):
    muonCounts = np.zeros(num)
    for i in range(num):
        proton = generateRandomPrimary(500000)
        muons = generateShower(proton)
        muonCounts[i] = len(muons)
    plt.hist(muonCounts,bins=muonCounts.max())
    plt.show()



if __name__ == '__main__':
    manyShowers(100)
