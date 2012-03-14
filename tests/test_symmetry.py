# -*- coding: utf-8 -*-
"""
test_symmetry
=============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import Symmetry,Shape,symmetry_generator2

def test_Symmetry_1D():
    h = Shape(1,np.array([2.0]))

    R = np.eye(1)
    t = np.zeros(1)
    s0 = Symmetry(1,"Bravais",h,R,t)
    is0 = s0.inverse()
    ms0 = s0 * is0

    R = -1.0 * np.eye(1)
    s1 = Symmetry(1,"Bravais",h,R,t)
    is1 = s1.inverse()
    ms1 = s1 * is1

    s2 = Symmetry(1,"Cartesian",h,R,t)
    is2 = s2.inverse()
    ms2 = s2 * is2

    print s0.__dict__
    print is0.__dict__
    print "np.eye(1) is expected:"
    print ms0.__dict__
    print s0 == is0

    print s1.__dict__
    print is1.__dict__
    print "np.eye(1) is expected:"
    print ms1.__dict__
    print s1 == is1

    print s2.__dict__
    print is2.__dict__
    print "np.eye(1) is expected:"
    print ms2.__dict__
    print s2 == is2

def test_Symmetry_2D():
    b = "Bravais"
    alpha = np.pi/3.0
    a1 = np.array([1.0,0.0])
    a2 = np.array([np.cos(alpha),np.sin(alpha)])
    h = Shape(2,a1,a2)

    R = np.eye(2)
    t = np.zeros(2)
    s0 = Symmetry(2,"Bravais",h,R,t)
    is0 = s0.inverse()
    ms0 = s0 * is0

    R = np.array([[0.0,-1.0],[1.0,-1.0]])
    t = np.zeros(2)
    s1 = Symmetry(2,"Bravais",h,R,t)
    is1 = s1.inverse()
    ms1 = s1 * is1

    R = np.eye(2)
    t = np.zeros(2)
    s2 = Symmetry(2,"Cartesian",h,R,t)
    is2 = s2.inverse()
    ms2 = s2 * is2

    R = np.array([[0.0,-1.0],[1.0,-1.0]])
    t = np.zeros(2)
    s3 = Symmetry(2,"Cartesian",h,R,t)
    is3 = s3.inverse()
    ms3 = s3 * is3

    print s0.__dict__
    print is0.__dict__
    print "np.eye(2) is expected:"
    print ms0.__dict__
    print s0 == is0

    print s1.__dict__
    print is1.__dict__
    print "np.eye(2) is expected:"
    print ms1.__dict__
    print s1 == is1

    print s2.__dict__
    print is2.__dict__
    print "np.eye(2) is expected:"
    print ms2.__dict__
    print s2 == is2

    print s3.__dict__
    print is3.__dict__
    print "np.eye(2) is expected:"
    print ms3.__dict__
    print s3 == is3

    symm = symmetry_generator2(17,b,h)
    for s in symm:
        print s.__dict__,"\n"
    if s1 in symm:
        print "s1 is in symm!"
    if is0 in symm:
        print "is0 is in symm!"
    if is1 in symm:
        print "is1 is in symm!"
    if is1 not in symm:
        print "is1 is not in symm!"

def test_Symmetry_3D():
    a1 = np.array([1.0,0.0,0.0])
    a2 = np.array([0.0,1.0,0.0])
    a3 = np.array([0.0,0.0,1.0])
    h = Shape(3,a1,a2,a3)

    R = np.eye(3)
    t = np.zeros(3)
    s0 = Symmetry(3,"Bravais",h,R,t)
    is0 = s0.inverse()
    ms0 = s0 * is0

    R = np.array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0]])
    t = np.array([0,0.5,0.5])
    s1 = Symmetry(3,"Bravais",h,R,t)
    is1 = s1.inverse()
    ms1 = s1 * is1

    R = np.eye(3)
    t = np.zeros(3)
    s2 = Symmetry(3,"Cartesian",h,R,t)
    is2 = s2.inverse()
    ms2 = s2 * is2

    R = np.array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0]])
    t = np.array([0,0.5,0.5])
    s3 = Symmetry(3,"Cartesian",h,R,t)
    is3 = s3.inverse()
    ms3 = s3 * is3

    print s0.__dict__
    print is0.__dict__
    print "np.eye(3) is expected:"
    print ms0.__dict__
    print s0 == is0

    print s1.__dict__
    print is1.__dict__
    print "np.eye(3) is expected:"
    print ms1.__dict__
    print s1 == is1

    print s2.__dict__
    print is2.__dict__
    print "np.eye(3) is expected:"
    print ms2.__dict__
    print s2 == is2

    print s3.__dict__
    print is3.__dict__
    print "np.eye(3) is expected:"
    print ms3.__dict__
    print s3 == is3

    symm = [s0,s1,s2,s3,ms0,ms1,ms2,ms3]
    if s1 in symm:
        print "s1 is in symm!"
    if is0 in symm:
        print "is0 is in symm!"
    if is1 not in symm:
        print "is1 is not in symm!"

def run_test():
    #test_Symmetry_1D()
    test_Symmetry_2D()
    #test_Symmetry_3D()

if __name__ == '__main__':
    run_test()
