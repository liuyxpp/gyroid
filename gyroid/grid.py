# -*- coding: utf-8 -*-
"""
gyroid.grid
===========

"""

import numpy as np
from numpy.linalg import inv

from .common import BRAVAIS,CARTESIAN
from .common import EPS,SMALL,LARGE

__all__ = ["Grid","wave_norm"]

class Grid(object):
    ''' A :class:`Grid` object represents a Discretized unit cell in reciprocal space.

    '''

    def __init__(self,ngrid,g,lowBZ=-1,highBZ=1):
        ''' Initiator of :class:`Grid`

        :param ngrid: contains number of grids in each dimension
        :type ngrid: `numpy.array` with integers
        :param g: a :class:`Group` object
        :param lowBZ: low bound to find G in BZ, default to -1
        :param highBZ: high bound to find G in BZ, default to 1

        '''

        self.dim = g.dim
        if np.size(ngrid) != g.dim:
            raise ValueError('Space dimension and the length of `ngrid` '
                             'not match.')
        self.N = ngrid
        self.shape = g.shape
        self.lowBZ = lowBZ
        self.highBZ = highBZ
        self.__create_waves(g)

    def to_BZ(self,G):
        ''' Shift a wave vector to the first Brillouin Zone (BZ).

        :param G: a wave vector
        :type G: `numpy.array`
        :returns: wave vector in BZ and its square magnitude.
        '''

        if np.size(G) != self.dim:
            raise ValueError('The wave vector and Dimension of '
                             'the grid not match in to_BZ.')

        if self.dim == 1:
            (i,) = G
            key = (i,)
            if not self.BZmap.has_key(key):
                ivec,ivec2 = self.__find_G_in_BZ(G)
                self.BZmap[key] = (ivec,ivec2)
            return self.BZmap[key]

        if self.dim == 2:
            (i,j) = G
            key = (i,j)
            if not self.BZmap.has_key(key):
                ivec,ivec2 = self.__find_G_in_BZ(G)
                self.BZmap[key] = (ivec,ivec2)
            return self.BZmap[key]

        if self.dim == 3:
            (i,j,k) = G
            key = (i,j,k)
            if not self.BZmap.has_key(key):
                ivec,ivec2 = self.__find_G_in_BZ(G)
                self.BZmap[key] = (ivec,ivec2)
            return self.BZmap[key]

    def is_wave_cancel(self,G,g):
        ''' Check whether the wave vector is canceled.

        A wave is canceled if and only if following conditions are met:

        1. Leaves **G** invariant (i.e. **G** . **R** == **G**)
        2. Produces a non-zero phase, such that **G** . **t** % 1.0 != 0

        :param G: a wave vector
        :type G: `numpy.array`
        :param g: a :class:`Group` object
        :rtype: True or False

        '''

        for i in np.arange(g.order):
            Gp = np.dot(G,g.symm[i].R)
            # Pseudo-Spetral method
            Gp,Gp2 = self.to_BZ(Gp)
            if np.all(np.abs(Gp-G) < EPS):
                phase = np.dot(G,g.symm[i].t) % 1.0
                # for cases phase=-1.0 - SMALL
                # (-1.0 - SMALL) % 1.0 ~ (1.0 - SMALL)
                if np.abs(phase-1.0) < EPS:
                    phase -= 1.0
                # wave canceled if phase not equal 0
                if np.abs(phase) > EPS:
                    return True
        return False

    @property
    def max_Gsq(self):
        ''' The maximun square of wave vector in the grid.'''

        length = SMALL

        if self.dim == 1:
            for (i,) in np.ndindex(self.N[0]):
                G = np.array([i])
                G,G2 = self.to_BZ(G)
                if G2 > length:
                    length = G2

        if self.dim == 2:
            for (i,j) in np.ndindex(self.N[0],self.N[1]):
                G = np.array([i,j])
                G,G2 = self.to_BZ(G)
                if G2 > length:
                    length = G2

        if self.dim == 3:
            for (i,j,k) in np.ndindex(self.N[0],self.N[1],self.N[2]):
                G = np.array([i,j,k])
                G,G2 = self.to_BZ(G)
                if G2 > length:
                    length = G2

        return length

    def __find_G_in_BZ(self,G):
        low = self.lowBZ
        high = self.highBZ
        G_try = G
        G_min = G
        Gsq_min = LARGE

        if self.dim == 1:
            for i in np.arange(high,low-1,-1):
                G_try = G + np.array([i]) * self.N
                Gsq = wave_norm(G_try,self.shape)
                if Gsq < Gsq_min:
                    Gsq_min, G_min = Gsq, G_try

        if self.dim == 2:
            for i in np.arange(high,low-1,-1):
                for j in np.arange(high,low-1,-1):
                    G_try = G + np.array([i,j]) * self.N
                    Gsq = wave_norm(G_try,self.shape)
                    if Gsq < Gsq_min:
                        Gsq_min, G_min = Gsq, G_try

        if self.dim == 3:
            for i in np.arange(high,low-1,-1):
                for j in np.arange(high,low-1,-1):
                    for k in np.arange(high,low-1,-1):
                        G_try = G + np.array([i,j,k]) * self.N
                        Gsq = wave_norm(G_try,self.shape)
                        if Gsq < Gsq_min:
                            Gsq_min, G_min = Gsq, G_try
        return G_min,Gsq_min

    def __create_waves(self,g):
        #Gsq_max = self.max_Gabs * self.max_Gabs
        #G_max = np.zeros(self.dim)
        #for i in np.arange(self.dim):
        #    aa = np.sqrt(np.dot(self.shape.h[i],self.shape.h[i]))
        #    G_max[i] = int(self.max_Gabs * aa / (2.0 * np.pi)) + 1

        # Calculate number of effective waves
        # Pseudo-Spectral method
        w,G2 = None,None
        self.BZmap = {}
        if self.dim == 1:
            for i in np.arange(self.N[0]):
                ivec,ivec2 = self.__find_G_in_BZ(np.array([i]))
                self.BZmap[(i,)] = (ivec,ivec2)
                if not self.is_wave_cancel(ivec,g):
                    if w is None:
                        w = np.array([ivec])
                        G2 = np.array([ivec2])
                    else:
                        w = np.append(w,[ivec],axis=0)
                        G2 = np.append(G2,[ivec2],axis=0)

        if self.dim == 2:
            for (i,j) in np.ndindex(self.N[0],self.N[1]):
                ivec,ivec2 = self.__find_G_in_BZ(np.array([i,j]))
                self.BZmap[(i,j)] = (ivec,ivec2)
                if not self.is_wave_cancel(ivec,g):
                    if w is None:
                        w = np.array([ivec])
                        G2 = np.array([ivec2])
                    else:
                        w = np.append(w,[ivec],axis=0)
                        G2 = np.append(G2,[ivec2],axis=0)

        if self.dim == 3:
            for (i,j,k) in np.ndindex(self.N[0],self.N[1],self.N[2]):
                ivec,ivec2 = self.__find_G_in_BZ(np.array([i,j,k]))
                self.BZmap[(i,j,k)] = (ivec,ivec2)
                if not self.is_wave_cancel(ivec,g):
                    if w is None:
                        w = np.array([ivec])
                        G2 = np.array([ivec2])
                    else:
                        w = np.append(w,[ivec],axis=0)
                        G2 = np.append(G2,[ivec2],axis=0)

        # Sort G2 in ascending order, returned the corresponding indices
        w = w.T
        self.Nw = np.size(w,1)
        ind = np.argsort(G2)
        self.waves = w[:,ind]
        self.Gsq = G2[ind]

def wave_norm(G,shape):
    ''' Calculate the square magnitude of a wave vector **G** ^2.

    :param G: a wave vector with 1, 2, or 3 elements
    :type G: 1D `numpy.array`
    :param shape: a :class:`Shape` object
    :return: **G** ^2
    :rtype: double

    '''

    v = np.dot(G,shape.g)
    return np.dot(v,v)

