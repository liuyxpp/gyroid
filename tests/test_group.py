# -*- coding: utf-8 -*-
"""
test_group
============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import Group,Shape

def test_Group1():
    b = "Bravais"
    h = Shape(1,np.array([2.0]))

    g1 = Group(1,b,h,"P1")
    #g2 = Group(1,"P2")
    g3 = Group(1,b,h,ITA_number = 2)
    #g4 = Group(1,ITA_number = 0)
    print g1.__dict__,"\n"
    for s in g1.symm:
        print s.__dict__,"\n"
    print g3.__dict__,"\n"
    for s in g3.symm:
        print s.__dict__,"\n"

def test_Group2():
    b = "Bravais"
    alpha = np.pi/3.0
    a1 = np.array([1.0,0.0])
    a2 = np.array([np.cos(alpha),np.sin(alpha)])
    h = Shape(2,a1,a2)

    g1 = Group(2,b,h,"P6mm")
    #g2 = Group(2,b,h,"P2")
    #g3 = Group(2,b,h,ITA_number = 17)
    #g4 = Group(2,b,h,ITA_number = 18)
    print g1.__dict__,"\n"
    for s in g1.symm:
        print s.__dict__,"\n"
    #print g3.__dict__,"\n"
    #for s in g3.symm:
    #    print s.__dict__,"\n"

def test_Group3():
    b = "Bravais"
    a1 = np.array([2.0,0.0,0.0])
    a2 = np.array([0.0,1.0,0.0])
    a3 = np.array([0.0,0.0,3.0])
    h = Shape(3,a1,a2,a3)

    #g1 = Group(3,b,h,"P6mm")
    #g2 = Group(3,b,h,"P2")
    g3 = Group(3,b,h,ITA_number = 230)
    #g4 = Group(3,b,h,ITA_number = 238)
    print g3.__dict__,"\n"
    for s in g3.symm:
        print s.__dict__,"\n"
    #print g3.__dict__,"\n"
    #for s in g3.symm:
    #    print s.__dict__,"\n"

def run_test():
    #test_Group1()
    #test_Group2()
    test_Group3()

if __name__ == '__main__':
    run_test()
