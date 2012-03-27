# -*- coding: utf-8 -*-
"""
test_star
=========

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
#import scipy.io
#import matplotlib.pyplot as plt
#from mayavi import mlab

from gyroid import UnitCell,Group,Grid,Basis
from gyroid.common import *
from gyroid import render_structure_1d
from gyroid import render_structure_2d
from gyroid import render_structure_3d

def test_Basis1():
    b = "Bravais"
    uc = UnitCell(1,LAMELLAR,np.array([4.0]))
    N1 = 128

    g1 = Group(1,b,uc.shape,"P1")
    gd1 = Grid(np.array([N1]),g1)
    bs = Basis(g1,gd1)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    c = 1.0
#    c = np.zeros(bs.N)
#    c[0] = 1.0 * N1
#    c[1] = 6.0
#    c[2] = 4.0
#    c[3] = 0.5
#    c[4] = 3.0
    render_structure_1d(bs,gd1,N1,c)

    g2 = Group(1,b,uc.shape,"P-1")
    gd2 = Grid(np.array([N1]),g2)
    bs2 = Basis(g2,gd2)
    print bs2.__dict__
    for s in bs2.stars:
        print s.__dict__
    c = 1.0
#    c = np.zeros(bs2.N)
#    c[0] = 1.0 * N1
#    c[1] = 6.0
#    c[2] = 1.0
#    c[3] = 3.0
#    c[4] = 0.2
    render_structure_1d(bs2,gd2,N1,c)

def test_Basis2():
    b = "Bravais"
    N1,N2 = 64,64

    uc = UnitCell(2,HEXAGONAL,np.array([1.0]))
    g1 = Group(2,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([64,64]),g1)
    bs = Basis(g1,gd1)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    c = np.zeros(bs.N)
    c[0] = 1.0*N1*N2
    c[1] = 18.0
    c[2] = 1.0
    render_structure_2d(bs,gd1,N1,N2,c)

def test_Basis3():
    b = "Bravais"

    N1,N2,N3 = 16,16,16
    uc = UnitCell(3,HEXAGONAL,np.array([1.0,1.0]))
    g = Group(3,b,uc.shape,"P6mm")
    gd = Grid(np.array([N1,N2,N3]),g)
    bs = Basis(g,gd)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    c = 1.0
#    c = np.zeros(bs.N)
#    c[0] = 1.0
#    c[1] = 1.0
#    c[2] = 1.0
#    c[3] = 1.0
#    c[4] = 1.0
    render_structure_3d(bs,gd,N1,N2,N3,c)

#    N1,N2,N3 = 32,32,32
#    uc = UnitCell(3)
#    g2 = Group(3,b,uc.shape,"Ia-3d")
#    gd2 = Grid(np.array([32,32,32]),g2)
#    bs2 = Basis(g2,gd2)
#    print bs2.__dict__
#    for s in bs2.stars:
#        print s.__dict__
#    c = np.zeros(bs2.N)
#    c[0] = 1
#    c[1] = 0.016
#    c[2] = 0.0042
#    c[3] = 0.00006
#    render_structure_3d(bs2,gd2,N1,N2,N3,c)

def run_test():
    test_Basis1()
    #test_Basis2()
    #test_Basis3()

if __name__ == '__main__':
    run_test()
