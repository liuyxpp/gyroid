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

    """
    A representation of a symmetry element in a group.

    The basis type of a symmetry element should be either 'Cartesian' or 'Bravais'.
    A symmetry element contains a point group matrix and a translational vector.

    """

    def __init__(self,dim,basis_type,shape_matrix,R=None,t=None):
        """
        *dim can be 1, 2 and 3
        *basis_type should be "Bravais" or "Cartesian"
        *shape_matrix is a shape matrix in Real space
        *R the point group symmetry matrix
        *t the translational symmetry vector. For __eq__ functioning correct, it must lie in the first unit cell.
        """
        self.dim = dim
        self.type = basis_type
        self.shape = shape_matrix
        if R is not None and np.size(R) == dim*dim:
            self.R = R
        else:
            raise ValueError("Dimensian and R not match")
        if t is not None and np.size(t) == dim:
            self.t = shape_matrix.shift(t,basis_type)
        else:
            raise ValueError("Dimensian and t not match")

    def inverse(self):
        iR = inv(self.R)
        it = -1.0 * np.dot(iR,self.t)
        it = self.shape.shift(it,self.type)
        return Symmetry(self.dim,self.type,self.shape,iR,it)

    def __mul__(self,rhs):
        if self.dim != rhs.dim or self.type != rhs.type:
            raise ValueError("Incompatible symmetries for multiplication.")
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


