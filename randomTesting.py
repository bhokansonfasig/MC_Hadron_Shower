"""File for testing distributions of randomly generated numbers"""
import numpy as np
import matplotlib.pyplot as plt
from interactions import rotate3D


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
    print(np.mean(a),np.mean(b),np.mean(c))
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
    print(np.mean(a),np.mean(b),np.mean(c))
    print(np.mean(s))
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
    a = []
    b = []
    c = []
    for _ in range(size):
        rands = np.random.random_sample(2)-.5
        if rands.sum()>.5 or rands.sum()<-.5:
            continue
        a.append(rands[0])
        b.append(rands[1])
        c.append(rands[0-rands.sum()])
    print(len(a)/size)
    print(np.mean(a),np.mean(b),np.mean(c))
    for x in [a,b,c]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test8():
    """Generation of three random vectors in a plane of fixed length,
    summing to zero"""
    size = 10000
    amag = 2
    bmag = 3
    cmag = 4
    ax = []
    ay = []
    bx = []
    by = []
    cx = []
    cy = []
    signs = []
    xsum = []
    ysum = []
    for i in range(size):
        xi = np.random.random_sample()*2*np.pi
        sign = np.random.random_sample()
        pm = int(sign<0.5)*2-1
        theta_a = xi
        theta_b = np.pi+xi-np.arccos((amag**2+bmag**2-cmag**2)/2/amag/bmag)
        theta_c = np.pi+xi+np.arccos((amag**2+cmag**2-bmag**2)/2/amag/cmag)
        ax.append(pm*amag*np.cos(theta_a))
        ay.append(pm*amag*np.sin(theta_a))
        bx.append(pm*bmag*np.cos(theta_b))
        by.append(pm*bmag*np.sin(theta_b))
        cx.append(pm*cmag*np.cos(theta_c))
        cy.append(pm*cmag*np.sin(theta_c))
        signs.append(pm)
        xsum = ax[i]+bx[i]+cx[i]
        ysum = ay[i]+by[i]+cy[i]
    for x in [ax,ay,bx,by,cx,cy,xsum,ysum,signs]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


def test9():
    """Generation of three random vectors in 3D of fixed length,
    summing to zero"""
    size = 10000
    amag = 2
    bmag = 3
    cmag = 4
    ax = []
    ay = []
    az = []
    bx = []
    by = []
    bz = []
    cx = []
    cy = []
    cz = []
    signs = []
    xsum = []
    ysum = []
    zsum = []
    for i in range(size):
        phi = np.random.random_sample()*2*np.pi
        theta = np.arccos(np.random.random_sample()*2-1)
        pm = int(np.random.random_sample()<0.5)*2-1
        xi_a = 0
        xi_b = np.pi-np.arccos((amag**2+bmag**2-cmag**2)/2/amag/bmag)
        xi_c = np.pi+np.arccos((amag**2+cmag**2-bmag**2)/2/amag/cmag)
        a = [0,0,amag]
        b = [bmag*np.sin(xi_b),0,bmag*np.cos(xi_b)]
        c = [cmag*np.sin(xi_c),0,cmag*np.cos(xi_c)]
        a2 = rotate3D(a,theta,phi)
        b2 = rotate3D(b,theta,phi)
        c2 = rotate3D(c,theta,phi)
        ax.append(pm*a2[0])
        ay.append(pm*a2[1])
        az.append(pm*a2[2])
        bx.append(pm*b2[0])
        by.append(pm*b2[1])
        bz.append(pm*b2[2])
        cx.append(pm*c2[0])
        cy.append(pm*c2[1])
        cz.append(pm*c2[2])
        signs.append(pm)
        xsum.append(ax[i]+bx[i]+cx[i])
        ysum.append(ay[i]+by[i]+cy[i])
        zsum.append(az[i]+bz[i]+cz[i])
    for x in [ax,ay,az,bx,by,bz,cx,cy,cz,xsum,ysum,zsum,signs]:
        plt.figure()
        plt.hist(x,bins=50)
        plt.show()


if __name__ == '__main__':
    # test0()
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # test6()
    # test7()
    # test8()
    test9()
