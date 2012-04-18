# -*- coding: utf-8 -*-
"""
gyroid.space_group
==================

"""

import numpy as np
from numpy.linalg import inv

from .symmetry import Symmetry

__all__ = ["symmetry_generator1",
           "symmetry_generator2",
           "symmetry_generator3"
          ]

def symmetry_generator1(ITA_number,b,h):
    ''' Output 1D space point group symmetry generators.

    There are **2** 1D space groups.
    Currently, all 1D space groups are supported::

       [1, 2]

    :param ITA_number: 1 for P1, 2 for P-1
    :type ITA_number: integer
    :param b: basis type
    :type b: string
    :param h: the shape of a unit cell
    :type h: :class:`Shape`
    :returns: the minimun number of symmetries which can be further expanded to the full set of point group symmetries
    :rtype: a list of :class:`Symmetry` objects

    '''

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
        raise ValueError('ITA number is not 1D.')

def symmetry_generator2(ITA_number,b,h):
    ''' Output 2D space point group symmetry generators.

    There are **17** 2D space groups.
    Currently, Only following 2D space groups are supported::

       [1,2,3,4,5,6,17]

    :param ITA_number: a sequential number as given in the International Tables for Crystallography, Vol. A, allowed range `[1,17]`
    :type ITA_number: integer
    :param b: basis type
    :type b: string
    :param h: the shape of a unit cell
    :type h: :class:`Shape`
    :returns: the minimun number of symmetries which can be further expande to the full set of point group symmetries
    :rtype: a list of :class:`Symmetry` objects

    '''

    I = np.eye(2)
    t0 = np.zeros(2)
    sI = Symmetry(2,b,h,I,t0)
    # q = 1/4, r = 3/4, s= 1/2, o = 1/8
    # Rotation
    R200 = np.array([[-1., 0.],
                     [ 0.,-1.]])

    R20y = np.array([[-1.,0.],
                     [ 0.,1.]])

    R300 = np.array([[0.,-1.],
                     [1.,-1.]])

    Rmxx = np.array([[ 0.,-1.],
                     [-1., 0.]])

    Rmx0 = np.array([[1., 0.],
                     [0.,-1.]])
    # Translation
    t0s = np.array([.0,.5])
    ts0 = np.array([.5,.0])
    tss = np.array([.5,.5])

    if ITA_number == 1:
        # P1, Oblique
        return [sI]

    if ITA_number == 2:
        # P2, Oblique
        return [sI,
                Symmetry(2,b,h,R200,t0)
               ]

    if ITA_number == 3:
        # Pm, Rectangular
        return [sI,
                Symmetry(2,b,h,R20y,t0)
               ]

    if ITA_number == 4:
        # Pg, Rectangular
        return [sI,
                Symmetry(2,b,h,R20y,t0s)
               ]

    if ITA_number == 5:
        # cm, Rectangular
        return [sI,
                Symmetry(2,b,h,R20y,t0),
                Symmetry(2,b,h,I,tss)
               ]

    if ITA_number == 6:
        # P2mm, Rectangular
        return [sI,
                Symmetry(2,b,h,R20y,t0),
                Symmetry(2,b,h,R200,t0)
               ]

    if ITA_number == 17:
        # "P6mm 2D", Hexagonal
        return [sI,
                Symmetry(2,b,h,R300,t0),
                Symmetry(2,b,h,R200,t0),
                Symmetry(2,b,h,Rmxx,t0)
               ]

    raise ValueError('ITA number not supported.')

def symmetry_generator3(ITA_number,b,h):
    ''' Output 3D space point group symmetry generators.

    There are **230** 3D space groups. The symmetry generators can
    be found at http://cryst.ehu.es/cryst/get_gen.html.
    Currently, the supported 3D space groups are::

       [183, 229, 230]

    :param ITA_number: a sequential number as given in the International Tables for Crystallography, Vol. A, allowed range `[1,17]`
    :type ITA_number: integer
    :param b: basis type
    :type b: string
    :param h: the shape of a unit cell
    :type h: :class:`Shape`
    :returns: the minimun number of symmetries which can be further expande to the full set of point group symmetries
    :rtype: a list of :class:`Symmetry`

    '''

    I = np.eye(3)
    t0 = np.zeros(3)
    sI = Symmetry(3,b,h,I,t0)
    # q = 1/4, r = 3/4, s= 1/2, o = 1/8
    # Rotation
    R300z = np.array([[0.,-1.,0.],
                      [1.,-1.,0.],
                      [0., 0.,1.]])

    R200z = np.array([[-1., 0.,0.],
                      [ 0.,-1.,0.],
                      [ 0., 0.,1.]])

    R20y0 = np.array([[-1.,0., 0.],
                      [ 0.,1., 0.],
                      [ 0.,0.,-1.]])

    Rmxxz = np.array([[ 0.,-1.,0.],
                      [-1., 0.,0.],
                      [ 0., 0.,1.]])

    R2q0z = np.array([[-1., 0.,0.],
                      [ 0.,-1.,0.],
                      [ 0., 0.,1.]])

    R20yq = np.array([[-1.,0., 0.],
                      [ 0.,1., 0.],
                      [ 0.,0.,-1.]])

    R3xxx = np.array([[0.,0.,1.],
                      [1.,0.,0.],
                      [0.,1.,0.]])

    R2xxo = np.array([[0.,1., 0.],
                      [1.,0., 0.],
                      [0.,0.,-1.]])

    R2xx0 = np.array([[0.,1., 0.],
                      [1.,0., 0.],
                      [0.,0.,-1.]])

    R1000 = np.array([[-1., 0., 0.],
                      [ 0.,-1., 0.],
                      [ 0., 0.,-1.]])
    # Translation
    Ts0s = np.array([.5,.0,.5])
    T0ss = np.array([.0,.5,.5])
    Trqq = np.array([.75,.25,.25])
    Tsss = np.array([.5,.5,.5])

    if ITA_number == 183:
        # "P6mm 3D", Hexagonal
        return [sI,
                Symmetry(3,b,h,R300z,t0),
                Symmetry(3,b,h,R200z,t0),
                Symmetry(3,b,h,Rmxxz,t0)
               ]
    if ITA_number == 229:
        # "Im-3m 3D", Cubic
        return [sI,
                Symmetry(3,b,h,R200z,t0),
                Symmetry(3,b,h,R20y0,t0),
                Symmetry(3,b,h,R3xxx,t0),
                Symmetry(3,b,h,R2xx0,t0),
                Symmetry(3,b,h,R1000,t0),
                Symmetry(3,b,h,I,Tsss)
               ]

    if ITA_number == 230:
        # "Ia-3d 3D", Cubic
        return [sI,
                Symmetry(3,b,h,R2q0z,Ts0s),
                Symmetry(3,b,h,R20yq,T0ss),
                Symmetry(3,b,h,R3xxx,t0),
                Symmetry(3,b,h,R2xxo,Trqq),
                Symmetry(3,b,h,R1000,t0),
                Symmetry(3,b,h,I,Tsss)
               ]
    else:
        raise ValueError('ITA number not supported.')

