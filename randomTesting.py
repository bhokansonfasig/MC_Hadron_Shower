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
        a[i] = np.sqrt(4*rands[0])
        b[i] = np.sqrt(4*(rands[1]-rands[0]))
        c[i] = np.sqrt(4*(1-rands[1]))
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


if __name__ == '__main__':
    # test0()
    # test1()
    # test2()
    test3()