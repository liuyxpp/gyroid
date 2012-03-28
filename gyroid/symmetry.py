# -*- coding: utf-8 -*-
"""
gyroid.symmetry
===============

"""

import numpy as np
from numpy.linalg import inv

from .common import EPS,BRAVAIS,CARTESIAN

__all__ = ["Symmetry"]

class Symmetry(object):
    ''' A :class:`Symmetry` object represents a point group symmetry element in a group.

    The basis type of a symmetry element should be either `BRAVAIS` or
    `CARTESIAN`.

    A symmetry element contains a point group matrix **R** and a translational
    vector **t**.

    '''

    def __init__(self,dim,basis_type,shape,R=None,t=None):
        '''
        :param dim: space dimension, can be 1, 2 or 3
        :type dim: integer
        :param basis_type: should be either `BRAVAIS` or `CARTESIAN`
        :type basis_type: string
        :param shape: the shape of a unit cell
        :type shape: :class:`Shape`
        :param R: the point group symmetry matrix
        :type R: `numpy.array`
        :param t: the translational symmetry vector. For __eq__ working
        correctly, it must lie in the first unit cell.
        :type t: 1D `numpy.array`

        '''

        self.dim = dim
        self.type = basis_type
        self.shape = shape
        if R is not None and np.size(R) == dim*dim:
            self.R = R
        else:
            raise ValueError('Dimensian and R not match')
        if t is not None and np.size(t) == dim:
            self.t = self.shape.shift(t,basis_type)
        else:
            raise ValueError('Dimensian and t not match')

    def inverse(self):
        ''' Inverse the symmetry element (both **R** and **t**)

        :return: the inverse symmetry element of the current object
        :rtype: :class:`Symmetry`

        '''
        iR = inv(self.R)
        it = -1.0 * np.dot(iR,self.t)
        it = self.shape.shift(it,self.type)
        return Symmetry(self.dim,self.type,self.shape,iR,it)

    def __mul__(self,rhs):
        if self.dim != rhs.dim or self.type != rhs.type:
            raise ValueError('Incompatible symmetries for multiplication.')
        R = np.dot(self.R,rhs.R)
        t = np.dot(self.R,rhs.t) + self.t
        t = self.shape.shift(t,self.type)
        return Symmetry(self.dim,self.type,self.shape,R,t)

    def __eq__(self,rhs):
        if self.dim != rhs.dim:
            return False
        if self.type != rhs.type:
            return False
        if np.any(np.abs(self.R - rhs.R) > EPS):
            return False
        if np.any(np.abs(self.t - rhs.t) > EPS):
            return False
        return True


