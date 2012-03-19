# -*- coding: utf-8 -*-
"""
test_group
============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import Group,Grid,Shape

def test_Grid1():
    b = "Bravais"
    h = Shape(1,np.array([2.0]))

    g1 = Group(1,b,h,"P1")
    gd1 = Grid(np.array([8]),g1)
    print gd1.__dict__,"\n"

    g3 = Group(1,b,h,ITA_number = 2)
    gd3 = Grid(np.array([9]),g3)
    print gd3.__dict__,"\n"

def test_Grid2():
    b = "Bravais"
    alpha = np.pi/3.0
    a1 = np.array([1.0,0.0])
    a2 = np.array([np.cos(alpha),np.sin(alpha)])
    h = Shape(2,a1,a2)

    g1 = Group(2,b,h,"P6mm")
    gd1 = Grid(np.array([4,4]),g1)
    print gd1.__dict__,"\n"

    gd3 = Grid(np.array([5,5]),g1)
    print gd3.__dict__,"\n"

    #g2 = Grid(2,b,h,"P2")
    #g3 = Grid(2,b,h,ITA_number = 17)
    #g4 = Grid(2,b,h,ITA_number = 18)
    #print g3.__dict__,"\n"
    #for s in g3.symm:
    #    print s.__dict__,"\n"

def test_Grid3():
    b = "Bravais"
    a1 = np.array([2.0,0.0,1.0])
    a2 = np.array([0.0,1.0,0.0])
    a3 = np.array([1.0,0.0,3.0])
    h = Shape(3,a1,a2,a3)

    g1 = Group(3,b,h,"P6mm")
    gd1 = Grid(np.array([3,3,3]),g1)
    print gd1.__dict__,"\n"

    gd2 = Grid(np.array([4,4,4]),g1)
    print gd2.__dict__,"\n"

    #g2 = Grid(3,b,h,"P2")
    #g3 = Grid(3,b,h,ITA_number = 230)
    #g4 = Grid(3,b,h,ITA_number = 238)
    #print g3.__dict__,"\n"
    #for s in g3.symm:
    #    print s.__dict__,"\n"

def run_test():
    #test_Grid1()
    #test_Grid2()
    test_Grid3()

if __name__ == '__main__':
    run_test()
