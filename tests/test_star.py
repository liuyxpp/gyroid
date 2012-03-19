# -*- coding: utf-8 -*-
"""
test_star
=========

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import UnitCell,Group,Grid,Star,Basis
from gyroid.common import *

def test_Star1():
    b = "Bravais"
    uc = UnitCell(1,LAMELLAR,np.array([2.0]))

    g1 = Group(1,b,uc.shape,"P1")
    gd1 = Grid(np.array([9]),g1)
    print gd1.Gsq
    s1 = Star(g1,gd1,gd1.Gsq[5])
    print s1.__dict__,"\n"

    g2 = Group(1,b,uc.shape,"P-1")
    gd2 = Grid(np.array([9]),g2)
    print gd2.Gsq
    s2 = Star(g1,gd1,gd1.Gsq[5])
    print s2.__dict__,"\n"

def test_Star2():
    b = "Bravais"
    uc = UnitCell(2,HEXAGONAL,np.array([2.0]))

    g1 = Group(2,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([128,128]),g1)
    print gd1.Gsq
    s1 = Star(g1,gd1,gd1.Gsq[7])
    print s1.__dict__

def test_Star3():
    b = "Bravais"

    uc = UnitCell(3,HEXAGONAL,np.array([2.0,5.0]))
    g1 = Group(3,b,uc.shape,"P6mm")
    gd1 = Grid(np.array([4,4,4]),g1)
    print gd1.Gsq
    s1 = Star(g1,gd1,gd1.Gsq[56])
    print s1.__dict__

#    uc = UnitCell(3)
#    g2 = Group(3,b,uc.shape,"Ia-3d")
#    gd2 = Grid(np.array([8,8,8]),g2)
#    print gd2.Gsq
#    s2 = Star(g2,gd2,gd2.Gsq[0])
#    print s2.__dict__
#    s3 = Star(g2,gd2,gd2.Gsq[49])
#    print s3.__dict__

def run_test():
    test_Star1()
    #test_Star2()
    #test_Star3()

if __name__ == '__main__':
    run_test()
