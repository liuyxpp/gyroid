# -*- coding: utf-8 -*-
"""
test_shape
============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import Shape

def test_Shape():
    twopi = 2.0*np.pi
    s1 = Shape(1,np.array([2.0]))
    print s1.l
    print s1.h
    print s1.g
    print s1.h*s1.g/twopi
    s2 = Shape(2,np.array([4.0,1.0]),np.array([2.0,3.0]))
    print s2.l
    print s2.h
    print s2.g
    print "The eye(2) matrix is expected."
    print np.dot(s2.g,s2.h.transpose())/twopi
    a1 = np.array([4.0,1.0,2.0])
    a2 = np.array([3.0,4.0,5.0])
    a3 = np.array([6.0,7.0,8.0])
    s3 = Shape(3,a1,a2,a3)
    print s3.l
    print s3.h
    print s3.g
    print "The eye(3) matrix is expected."
    print np.dot(s3.g,s3.h.transpose())/twopi


def run_test():
    # test_inverse() passed.
    test_Shape()

if __name__ == '__main__':
    run_test()
