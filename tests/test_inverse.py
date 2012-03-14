# -*- coding: utf-8 -*-
"""
test_inverse
============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from numpy.linalg import inv

def inverse(m):
    """ 1x1, 2x2, and 3x3 matrix inversion. """
    if np.size(m) == 1:
        return 1.0/m
    return inv(m)

def test_inverse():
    print "TEST inverse() for 3x3, 2x2 matrix and scalar..."
    a = np.array([[3,1,5],[1,0,8],[2,1,4]])
    ia = inverse(a)
    print "The resulted matrix should be eye(3) within machine precision."
    print np.dot(a,ia)
    b = np.array([[4,1],[2,3]])
    ib = inverse(b)
    print "The resulted matrix should be eye(2) within machine precision."
    print np.dot(b,ib)
    c = 2.0
    ic = inverse(c)
    print "The resulted value should be 1.0 within machine precision."
    print c * ic

def run_test():
    test_inverse()

if __name__ == '__main__':
    run_test()
