# -*- coding: utf-8 -*-
"""
gyroid.unitcell
===============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from numpy.linalg import inv

from .common import BRAVAIS,CARTESIAN

__all__ = ["Shape"]

class Shape(object):

    """
    Shape matrix constructed from unit vectors in Cartesian Coordinate.

    The Morse convention is used. That is each row in the shape
    matrix represents a unit vector, e.g. h = (a1,a2,a3) where a_i =
    (x_i,y_i,z_i) is the unit vector of the Bravis lattice in Cartesian
    Coordinate.
    """

    def __init__(self,dim,a1,a2=None,a3=None):
        self.dim = dim
        if dim==3 and np.size(a1)==3 and np.size(a2)==3 and np.size(a3)==3:
            self.m = np.array([a1,a2,a3])
        elif dim==2 and np.size(a1)==2 and np.size(a2)==2:
            self.m = np.array([a1,a2])
        elif dim==1 and np.size(a1)==1:
            self.m = np.array([a1])
        else:
            raise ValueError(
                "Dimension and the number of unit vector not match")

    def shift(self,t,basis_type):
        """
        Shift position or translation vector so as to lie in the first unit
        cell.
        """
        if basis_type == BRAVAIS:
            return t % 1.0
        elif basis_type == CARTESIAN:
            cc = np.dot(self.g,t) / (2.0 * np.pi)
            return np.dot(self.h, cc % 1.0)
        else:
            raise ValueError("Unrecognized basis type when shift a vector.")

    @property
    def h(self):
        """ The shape matrix in Real Space. """
        return self.m

    @property
    def g(self):
        """ The shape matrix in Reciprocal Space. """
        return 2.0 * np.pi * inv(self.m).T


