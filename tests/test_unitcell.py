# -*- coding: utf-8 -*-
"""
test_unitcell
=============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from gyroid import UnitCell
from gyroid.common import *

def test_UnitCell1():
    uc = UnitCell(1,LAMELLAR,np.array([2.0]))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(1,cell_param=np.array([3.0]))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(1)
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

def test_UnitCell2():
    uc = UnitCell(2,HEXAGONAL,np.array([2.0]))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(2,OBLIQUE,cell_param=np.array([3.0,4.0,np.pi/4.0]))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(2)
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

def test_UnitCell3():
    uc = UnitCell(3,HEXAGONAL,np.array([2.0,5.0]))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(3,TRICLINIC,
                  cell_param=np.array(
                      [3.0,4.0,5.0,np.pi/3.0,np.pi/4.0,np.pi/5.0]
                  ))
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

    uc = UnitCell(3)
    print uc.shape.__dict__
    print uc.shape.h
    print uc.shape.g

def run_test():
    #test_UnitCell1()
    #test_UnitCell2()
    test_UnitCell3()

if __name__ == '__main__':
    run_test()
