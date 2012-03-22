# -*- coding: utf-8 -*-
"""
test_star
=========

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mayavi import mlab

from gyroid import UnitCell,Group,Grid,Star,Basis
from gyroid.common import *

def test_Basis1():
    b = "Bravais"
    uc = UnitCell(1,LAMELLAR,np.array([2.0]))

    g1 = Group(1,b,uc.shape,"P1")
    gd1 = Grid(np.array([4]),g1)
    bs = Basis(g1,gd1)
    struct = bs.generate_structure(128,1.0)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    print struct
    plt.plot(struct)
    plt.show()

    g2 = Group(1,b,uc.shape,"P-1")
    gd2 = Grid(np.array([128]),g2)
    bs2 = Basis(g2,gd2)
    struct = bs2.generate_structure(128,1.0)
    print bs2.__dict__
    for s in bs2.stars:
        print s.__dict__
    print struct
    plt.plot(struct)
    plt.show()

def test_Basis2():
    b = "Bravais"

    uc = UnitCell(2,HEXAGONAL,np.array([2.0]))
    g1 = Group(2,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([4,4]),g1)
    bs = Basis(g1,gd1)
    struct = bs.generate_structure((128,128),1.0)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    plt.matshow(struct)
    plt.show()

#    uc = UnitCell(2) #,RECTANGULAR,np.array([2.0,np.sqrt(3.0)])
#    g2 = Group(2,b,uc.shape,"P6mm")
#    gd2 = Grid(np.array([8,8]),g2)
#    bs2 = Basis(g2,gd2)
#    struct = bs2.generate_structure((128,128),1.0)
#    print bs2.__dict__
#    for s in bs2.stars:
#        print s.__dict__
#    plt.matshow(struct)
#    plt.show()

def test_Basis3():
    b = "Bravais"
    N1,N2,N3 = 16,16,32

    uc = UnitCell(3,HEXAGONAL,np.array([1.0,2.0]))
    g1 = Group(3,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([4,4,4]),g1)
    bs = Basis(g1,gd1)
    struct = bs.generate_structure((N1,N2,N3),1.0)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__

    rx = np.zeros((N1,N2,N3))
    ry = np.zeros((N1,N2,N3))
    rz = np.zeros((N1,N2,N3))
    for (i,j,k) in np.ndindex(N1,N2,N3):
        x = (1.0*np.array([i,j,k])) / (N1-1,N2-1,N3-1)
        r = np.dot(x,uc.shape.h)
        rx[i,j,k],ry[i,j,k],rz[i,j,k] = r[0],r[1],r[2]
    print uc.shape.h

    scipy.io.savemat("struct3d.mat",{"x":rx,"y":ry,"z":rz,"s":struct})
    mlab.contour3d(rx,ry,rz,struct)
    img = mlab.screenshot()
    plt.imshow(img)
    plt.show()

#    uc = UnitCell(3)
#    g2 = Group(3,b,uc.shape,"Ia-3d")
#    gd2 = Grid(np.array([8,8,8]),g2)
#    print gd2.Gsq
#    s2 = Star(g2,gd2,gd2.Gsq[0])
#    print s2.__dict__
#    s3 = Star(g2,gd2,gd2.Gsq[49])
#    print s3.__dict__

def run_test():
    #test_Basis1()
    #test_Basis2()
    test_Basis3()

if __name__ == '__main__':
    run_test()
