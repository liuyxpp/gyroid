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

    g1 = Group(1,b,uc.shape,"P1")
    gd1 = Grid(np.array([4]),g1)
    bs = Basis(g1,gd1)
    render_structure_1d(bs,128,1.0)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__

    g2 = Group(1,b,uc.shape,"P-1")
    gd2 = Grid(np.array([128]),g2)
    bs2 = Basis(g2,gd2)
    render_structure_1d(bs2,128,1.0)
    print bs2.__dict__
    for s in bs2.stars:
        print s.__dict__

def test_Basis2():
    b = "Bravais"
    N1,N2 = 64,64

    uc = UnitCell(2,HEXAGONAL,np.array([2.0]))
    g1 = Group(2,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([32,32]),g1)
    bs = Basis(g1,gd1)
    print bs.__dict__
    for s in bs.stars:
        print s.__dict__
    render_structure_2d(bs,N1,N2,1.0)

def test_Basis3():
    b = "Bravais"

#    N1,N2,N3 = 16,16,32
#    uc = UnitCell(3,HEXAGONAL,np.array([1.0,2.0]))
#    g1 = Group(3,b,uc.shape,"P6mm")
#    gd1 = Grid(np.array([4,4,4]),g1)
#    bs = Basis(g1,gd1)
#    render_structure_3d(bs,N1,N2,N3,1.0)
#    print bs.__dict__
#    for s in bs.stars:
#        print s.__dict__

    N1,N2,N3 = 32,32,32
    uc = UnitCell(3)
    g2 = Group(3,b,uc.shape,"Ia-3d")
    gd2 = Grid(np.array([5,5,5]),g2)
    bs2 = Basis(g2,gd2)
    print bs2.__dict__
    for s in bs2.stars:
        print s.__dict__
    render_structure_3d(bs2,N1,N2,N3,1.0)

def run_test():
    #test_Basis1()
    test_Basis2()
    #test_Basis3()

if __name__ == '__main__':
    run_test()
