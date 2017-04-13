"""File for testing distributions of randomly generated numbers"""
import numpy as np
import matplotlib.pyplot as plt


def test0():
    """Incorrect generation of three random numbers with a fixed sum"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(2)
        a[i] = rands[0]
        b[i] = rands[1]
        c[i] = 3-a[i]-b[i]
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test1():
    """Incorrect generation of three random numbers with a fixed sum"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(3)
        a[i] = 4*rands[0]/rands.sum()
        b[i] = 4*rands[1]/rands.sum()
        c[i] = 4*rands[2]/rands.sum()
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test2():
    """Correct generation of three random numbers with a fixed sum"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(2)
        rands.sort()
        a[i] = 4*rands[0]
        b[i] = 4*(rands[1]-rands[0])
        c[i] = 4*(1-rands[1])
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test3():
    """Potentially bad way of generating vector components from magnitude"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(2)
        rands.sort()
        negativeSeeds = np.random.random_sample(3)
        signs = [1,1,1]
        for j in range(3):
            if negativeSeeds[j]<.5:
                signs[j] = -1
        a[i] = signs[0]*np.sqrt(4*rands[0])
        b[i] = signs[1]*np.sqrt(4*(rands[1]-rands[0]))
        c[i] = signs[2]*np.sqrt(4*(1-rands[1]))
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test4():
    """Desired distribution for vector components"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    cosThetas = np.random.random_sample(size)*2-1
    thetas = np.arccos(cosThetas)
    phis = np.random.random_sample(size)*2*np.pi
    for i in range(size):
        theta = thetas[i]
        phi = phis[i]
        a[i] = 2*np.sin(theta)*np.cos(phi)
        b[i] = 2*np.sin(theta)*np.sin(phi)
        c[i] = 2*np.cos(theta)
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test5():
    """Bad way of distributing vector components among three vectors"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    cosThetas = np.random.random_sample(2*size)*2-1
    thetas = np.arccos(cosThetas)
    phis = np.random.random_sample(2*size)*2*np.pi
    for i in range(size):
        theta1 = thetas[i]
        phi1 = phis[i]
        theta2 = thetas[size+i]
        phi2 = phis[size+i]
        a[i] = 0-(np.sin(theta1)*np.cos(phi1) + np.sin(theta2)*np.cos(phi2))
        b[i] = 0-(np.sin(theta1)*np.sin(phi1) + np.sin(theta2)*np.sin(phi2))
        c[i] = 0-(np.cos(theta1) + np.cos(theta2))
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test6():
    """Generation of three random numbers with a fixed sum of zero"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    s = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(2)
        rands.sort()
        a[i] = rands[0] - 1/3
        b[i] = rands[1]-rands[0] - 1/3
        c[i] = 1-rands[1] - 1/3
        s[i] = a[i]+b[i]+c[i]
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()
    plt.figure()
    plt.plot(s)
    plt.show()


def test7():
    """Generation of three random numbers with a fixed sum of zero"""
    size = 10000
    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    s = np.zeros(size)
    for i in range(size):
        rands = np.random.random_sample(2)
        rands.sort()
        rands -= .5
        a[i] = rands[0]+.5
        b[i] = rands[1]-rands[0]
        c[i] = .5-rands[1]
        s[i] = a[i]+b[i]+c[i]
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()
    plt.figure()
    plt.plot(s)
    plt.show()


if __name__ == '__main__':
    # test0()
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # test6()
    test7()
