# -*- coding: utf-8 -*-
"""
gyroid.unitcell
===============

Define a standard unit cell and its (real space) lattice basis vectors according to the cystal system.

"""

import numpy as np
from numpy.linalg import inv

from .common import BRAVAIS,CARTESIAN
from .common import CRYSTAL_SYSTEM1,CRYSTAL_SYSTEM2,CRYSTAL_SYSTEM3
from .common import LAMELLAR
from .common import SQUARE,RECTANGULAR,HEXAGONAL,OBLIQUE
from .common import CUBIC,TETRAGONAL,ORTHORHOMBIC,TRIGONAL,MONOCLINIC
from .common import TRICLINIC,DEFAULT

__all__ = ["UnitCell","Shape"]

class UnitCell(object):
    ''' A :class:`UnitCell` object contains the shape information of a unit
    cell.

    All standard crystal systems in 1D, 2D and 3D spaces are supported.

    To construct a unit cell of a standard crystal system, one should provide
    the space dimension, the name of the standard crystal system, and the cell
    parameters (the lengths of each unit vectors in the unit cell, and the angles
    between the unit vectors).


    The name of the standard crystal system is best given by first import the
    gyroid.common module, and use the constant therein. However, you can also
    input a string, such as 'Cubic' for 3D cubic crystal systems. If the name of
    the crystal system is not given, it will use the default crystal system for
    each space dimension, and the cell parameters given by the user will be
    discarded.

    The default unit cell in 1D space is `LAMELLAR`, unit vector length is 1.0.
    The default unit cell in 2D space is `SQUARE`, unit vector length is 1.0.
    The default unit cell in 3D space is `CUBIC`, unit vector length is also
    1.0.

    '''

    def __init__(self,dim,cell_type=DEFAULT,cell_param=None):
        ''' If cell_type is missing, cell_param has no effect!'''

        self.dim = dim
        if dim == 1:
            if cell_type in CRYSTAL_SYSTEM1:
                self.type = cell_type
                self.__standard_cell_1D(cell_type,cell_param)
            else:
                raise ValueError('Unkonwn crystal system for 1D space.')

        if dim == 2:
            if cell_type in CRYSTAL_SYSTEM2:
                self.type = cell_type
                self.__standard_cell_2D(cell_type,cell_param)
            else:
                raise ValueError('Unkonwn crystal system for 2D space.')

        if dim == 3:
            if cell_type in CRYSTAL_SYSTEM3:
                self.type = cell_type
                self.__standard_cell_3D(cell_type,cell_param)
            else:
                raise ValueError('Unkonwn crystal system for 3D space.')

        self.shape = self.__create_shape()

    def __create_shape(self):
        '''
        :var a: length of Bravais unit vector **a**
        :var b: length of Bravais unit vector **b**
        :var c: length of Bravais unit vector **c**
        :var alpha: angle between vectors **b** and **c**
        :var beta: angle between vectors **c** and **a**
        :var gamma: angle between vectors **a** and **b**
        :retrun: the shape of the unit cell
        :rtype: :class:`Shape`

        '''

        if self.dim == 1:
            return Shape(1,np.array([self.a]))

        if self.dim == 2:
            a = np.array([self.a, 0.0])
            b = np.array([self.b*np.cos(self.gamma),
                          self.b*np.sin(self.gamma)])
            return Shape(2,a,b)

        if self.dim == 3:
            a = np.array([self.a, 0.0, 0.0])
            b = np.array([self.b*np.cos(self.gamma),
                          self.b*np.sin(self.gamma),
                          0.0])
            cx = self.c*np.cos(self.beta)
            cy = self.c*(np.cos(self.alpha)-
                         np.cos(self.beta)*np.cos(self.gamma))
            cz = np.sqrt(self.c*self.c - cx*cx - cy*cy)
            c = np.array([cx,cy,cz])
            return Shape(3,a,b,c)

    def __standard_cell_1D(self,cell_type,cp):
        if cell_type == LAMELLAR:
            if np.size(cp) < 1:
                raise ValueError('Lamellar crystal requires 1 parameters.')
            self.a = cp[0]
        elif cell_type == DEFAULT:
            self.a = 1.0
        else:
            raise ValueError('Unknow 1D crystal system.')

    def __standard_cell_2D(self,cell_type,cp):
        pi2 = np.pi / 2.0
        pi3 = 2.0 * np.pi / 3.0
        if cell_type == SQUARE:
            if np.size(cp) < 1:
                raise ValueError('Square crystal requires 1 parameters.')
            self.a, self.b = cp[0], cp[0]
            self.gamma = pi2
        elif cell_type == RECTANGULAR:
            if np.size(cp) < 2:
                raise ValueError('Rectangular crystal '
                                 'requires 2 parameters.')
            self.a, self.b = cp[0], cp[1]
            self.gamma = pi2
        elif cell_type == HEXAGONAL:
            if np.size(cp) < 1:
                raise ValueError('Hexagonal crystal '
                                 'requires 1 parameters.')
            self.a, self.b = cp[0], cp[0]
            self.gamma = pi3
        elif cell_type == OBLIQUE:
            if np.size(cp) < 3:
                raise ValueError('Oblique crystal '
                                 'requires 3 parameters.')
            self.a, self.b = cp[0], cp[1]
            self.gamma = cp[2] * (np.pi / 180.0)
        elif cell_type == DEFAULT:
            self.a, self.b = 1.0, 1.0
            self.gamma = pi2
        else:
            raise ValueError('Unknow 2D crystal system.')

    def __standard_cell_3D(self,cell_type,cp):
        pi2 = np.pi / 2.0
        pi3 = 2.0 * np.pi / 3.0
        if cell_type == CUBIC:
            if np.size(cp) < 1:
                raise ValueError('Cubic crystal requires 1 parameters.')
            self.a, self.b, self.c = cp[0], cp[0], cp[0]
            self.alpha, self.beta, self.gamma = pi2, pi2, pi2
        elif cell_type == TETRAGONAL:
            if np.size(cp) < 2:
                raise ValueError('Tetragonal crystal '
                                 'requires 2 parameters.')
            self.a, self.b, self.c = cp[0], cp[0], cp[1]
            self.alpha, self.beta, self.gamma = pi2, pi2, pi2
        elif cell_type == ORTHORHOMBIC:
            if np.size(cp) < 3:
                raise ValueError('Orthorhombic crystal '
                                 'requires 3 parameters.')
            self.a, self.b, self.c = cp[0], cp[1], cp[2]
            self.alpha, self.beta, self.gamma = pi2, pi2, pi2
        elif cell_type == HEXAGONAL:
            if np.size(cp) < 2:
                raise ValueError('Hexagonal crystal '
                                 'requires 3 parameters.')
            self.a, self.b, self.c = cp[0], cp[0], cp[1]
            self.alpha, self.beta, self.gamma = pi2, pi2, pi3
        elif cell_type == TRIGONAL:
            if np.size(cp) < 2:
                raise ValueError('Trigonal crystal '
                                 'requires 3 parameters.')
            self.a, self.b, self.c = cp[0], cp[0], cp[0]
            radian = cp[1] * (np.pi / 180.0)
            self.alpha, self.beta, self.gamma = radian, radian, radian
        elif cell_type == MONOCLINIC:
            if np.size(cp) < 4:
                raise ValueError('Monoclinic crystal '
                                 'requires 3 parameters.')
            self.a, self.b, self.c = cp[0], cp[1], cp[2]
            radian = cp[3] * (np.pi / 180.0)
            self.alpha, self.beta, self.gamma = pi2, radian, pi2
        elif cell_type == TRICLINIC:
            if np.size(cp) < 6:
                raise ValueError('Triclinic crystal '
                                 'requires 6 parameters.')
            self.a, self.b, self.c = cp[0], cp[1], cp[2]
            radian1 = cp[3] * (np.pi / 180.0)
            radian2 = cp[4] * (np.pi / 180.0)
            radian3 = cp[5] * (np.pi / 180.0)
            self.alpha, self.beta, self.gamma = radian1, radian2, radian3
        elif cell_type == DEFAULT:
            self.a, self.b, self.c = 1.0, 1.0, 1.0
            self.alpha, self.beta, self.gamma = pi2, pi2, pi2
        else:
            raise ValueError('Unknown 3D crystal system.')


class Shape(object):
    ''' A :class:`Shape` object fully describes the shape of a unit cell.

    The heart of a :class:`Shape` object is a shape matrix which can be
    constructed from unit vectors of a unit cell in Cartesian Coordinate.

    The Morse convention is used to construct the shape matrix. That is each row
    in the shape matrix represents a unit vector, e.g.::

       h = [a1 a2 a3
            b1 b2 b3
            c1 c2 c3]

    where `a_i`, `b_i`, `c_i` are the *i* th elements of unit vectors **a**,
    **b**, and **c** in the Bravis lattice in real space.

    '''

    def __init__(self,dim,a1,a2=None,a3=None):
        self.dim = dim
        if dim==3 and np.size(a1)==3 and np.size(a2)==3 and np.size(a3)==3:
            self.m = np.array([a1,a2,a3])
        elif dim==2 and np.size(a1)==2 and np.size(a2)==2:
            self.m = np.array([a1,a2])
        elif dim==1 and np.size(a1)==1:
            self.m = np.array([a1])
        else:
            raise ValueError('Dimension and the number'
                             'of unit vector not match')
        self.gg = 2.0 * np.pi * inv(self.m).T
        self.ll = np.sqrt(np.sum(self.m**2,axis=1))

    def shift(self,t,basis_type):
        ''' Shift position or translation vector so as to lie in the first unit
        cell.

        :param t: translational or positional vector to be shifted.
        :type t: 1D `numpy.array`
        :param basis_type: the type of the unit vector, can be either `BRAVAIS` or `CARTESIAN`
        :type basis_type: string
        :return: shifted translational or positional vector
        :rtype: 1D `numpy.array`

        '''

        if basis_type == BRAVAIS:
            return t % 1.0
        elif basis_type == CARTESIAN:
            cc = np.dot(t,self.g) / (2.0 * np.pi)
            return np.dot(self.h, cc % 1.0)
        else:
            raise ValueError('Unrecognized basis type when shift a vector.')

    @property
    def l(self):
        ''' (a,b,c), the length vectors of each dimension. '''
        return self.ll

    @property
    def h(self):
        ''' The shape matrix in real space. '''

        return self.m

    @property
    def g(self):
        ''' The shape matrix in reciprocal space. '''

        return self.gg


