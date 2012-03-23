# -*- coding: utf-8 -*-
"""
gyroid.space_group
===============

"""

import numpy as np
from numpy.linalg import inv

from .symmetry import Symmetry

__all__ = ["symmetry_generator1",
           "symmetry_generator2",
           "symmetry_generator3"
          ]

def symmetry_generator1(ITA_number,b,h):
    I = np.eye(1)
    t0 = np.zeros(1)
    sI = Symmetry(1,b,h,I,t0)
    if ITA_number == 1:
        return [sI]
    elif ITA_number == 2:
        return [sI,
                Symmetry(1,b,h,-1.0*I,t0)
               ]
    else:
        raise ValueError("ITA number is not 1D.")

def symmetry_generator2(ITA_number,b,h):
    """
    There are 17 2D space groups.
    Currently, Only following 2D space groups are supported:
        [17]

    Variable
    --------
        ITA_number  -   a sequential number as given in the 
                        International Tables for Crystallography, 
                        Vol. A
        b           -   Basis type
        h           -   Shape instance that describes the unit cell

    Output
    ------
        A list of Symmetry instances that contains the minimun number of
        symmetries which can be further expanded to the full set of point 
        group symmetries.
    """
    I = np.eye(2)
    t0 = np.zeros(2)
    sI = Symmetry(2,b,h,I,t0)
    R300 = np.array([[0.0,-1.0],[1.0,-1.0]])
    R200 = np.array([[-1.0,0.0],[0.0,-1.0]])
    Rmxx = np.array([[0.0,-1.0],[-1.0,0.0]])
    if ITA_number == 17:
        # "P6mm 2D"
        return [sI,
                Symmetry(2,b,h,R300,t0),
                Symmetry(2,b,h,R200,t0),
                Symmetry(2,b,h,Rmxx,t0)
               ]
    else:
        raise ValueError("ITA number not supported.")

def symmetry_generator3(ITA_number,b,h):
    I = np.eye(3)
    t0 = np.zeros(3)
    sI = Symmetry(3,b,h,I,t0)
    # q = 1/4, r = 3/4, s= 1/2, o = 1/8
    # Rotation
    R300z = np.array([[0.0,-1.0,0.0],[1.0,-1.0,0.0],[0.0,0.0,1.0]])
    R200z = np.array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0]])
    Rmxxz = np.array([[0.0,-1.0,0.0],[-1.0,0.0,0.0],[0.0,0.0,1.0]])
    R2q0z = np.array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0]])
    R20yq = np.array([[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0]])
    R3xxx = np.array([[0.0,0.0,1.0],[1.0,0.0,0.0],[0.0,1.0,0.0]])
    R2xxo = np.array([[0.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,-1.0]])
    R1000 = np.array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,-1.0]])
    # Translation
    Ts0s = np.array([0.5,0.0,0.5])
    T0ss = np.array([0.0,0.5,0.5])
    Trqq = np.array([0.75,0.25,0.25])
    Tsss = np.array([0.5,0.5,0.5])
    if ITA_number == 183:
        # "P6mm 3D"
        return [sI,
                Symmetry(3,b,h,R300z,t0),
                Symmetry(3,b,h,R200z,t0),
                Symmetry(3,b,h,Rmxxz,t0)
               ]
    if ITA_number == 230:
        # "Ia-3d 3D"
        return [sI,
                Symmetry(3,b,h,R2q0z,Ts0s),
                Symmetry(3,b,h,R20yq,T0ss),
                Symmetry(3,b,h,R3xxx,t0),
                Symmetry(3,b,h,R2xxo,Trqq),
                Symmetry(3,b,h,R1000,t0),
                Symmetry(3,b,h,I,Tsss)
               ]
    else:
        raise ValueError("ITA number not supported.")

