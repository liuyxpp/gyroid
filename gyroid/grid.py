# -*- coding: utf-8 -*-
"""
gyroid.grid
===============

"""

import numpy as np
from numpy.linalg import inv

from .common import BRAVAIS,CARTESIAN
from .common import EPS,SMALL,LARGE

__all__ = ["Grid","wave_norm"]

class Grid(object):
    """
    Discrete form of a unit cell in reciprocal space.
    """

    def __init__(self,ngrid,g):
        """
        ngrid must be a length dim vector containing positive integers.
        ngrid must be a numpy.ndarray.
        g is a Group instance.
        """
        self.dim = g.dim
        if np.size(ngrid) != g.dim:
            raise ValueError("Dimension and ngrid not match.")
        self.N = ngrid
        self.shape = g.shape
        self.__create_waves(g)

    def to_BZ(self,G):
        if np.size(G) != self.dim:
            raise ValueError("""The wave vector and Dimension of
                             the grid not match in to_BZ.""")
        low = -1
        high = 1
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
        return G_min

    def is_wave_cancel(self,G,g):
        """
        A wave is canceled if and only if following conditions are met:
            1) Leaves G invariant (i.e. G.R == G), and
            2) Produces a non-zero phase, such that G.t % 1.0 != 0
        """
        for i in np.arange(g.order):
            Gp = np.dot(G,g.symm[i].R)
            # Pseudo-Spetral method
            Gp = self.to_BZ(Gp)
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
    def max_Gabs(self):
        length = SMALL

        if self.dim == 1:
            for (i,) in np.ndindex(self.N[0]):
                G = np.array([i])
                tmp = wave_norm(self.to_BZ(G),self.shape)
                if tmp > length:
                    length = tmp

        if self.dim == 2:
            for (i,j) in np.ndindex(self.N[0],self.N[1]):
                G = np.array([i,j])
                tmp = wave_norm(self.to_BZ(G),self.shape)
                if tmp > length:
                    length = tmp

        if self.dim == 3:
            for (i,j,k) in np.ndindex(self.N[0],self.N[1],self.N[2]):
                G = np.array([i,j,k])
                tmp = wave_norm(self.to_BZ(G),self.shape)
                if tmp > length:
                    length = tmp

        return np.sqrt(length)

    def __create_waves(self,g):
        #Gsq_max = self.max_Gabs * self.max_Gabs
        #G_max = np.zeros(self.dim)
        #for i in np.arange(self.dim):
        #    aa = np.sqrt(np.dot(self.shape.h[i],self.shape.h[i]))
        #    G_max[i] = int(self.max_Gabs * aa / (2.0 * np.pi)) + 1

        # Calculate number of effective waves
        # Pseudo-Spectral method
        if self.dim == 1:
            n = 0
            for (i,) in np.ndindex(self.N[0]):
                ivec = self.to_BZ(np.array([i]))
                if not self.is_wave_cancel(ivec,g):
                    n += 1
            Nw = n
            w = np.zeros((self.dim,Nw))
            G2 = np.zeros(Nw)
            n = 0
            for (i,) in np.ndindex(self.N[0]):
                ivec = self.to_BZ(np.array([i]))
                if not self.is_wave_cancel(ivec,g):
                    w[:,n] = ivec
                    G2[n] = wave_norm(ivec,self.shape)
                    n += 1

        if self.dim == 2:
            n = 0
            for (i,j) in np.ndindex(self.N[0],self.N[1]):
                ivec = self.to_BZ(np.array([i,j]))
                if not self.is_wave_cancel(ivec,g):
                    n += 1
            Nw = n
            w = np.zeros((self.dim,Nw))
            G2 = np.zeros(Nw)
            n = 0
            for (i,j) in np.ndindex(self.N[0],self.N[1]):
                ivec = self.to_BZ(np.array([i,j]))
                if not self.is_wave_cancel(ivec,g):
                    w[:,n] = ivec
                    G2[n] = wave_norm(ivec,self.shape)
                    n += 1

        if self.dim == 3:
            n = 0
            for (i,j,k) in np.ndindex(self.N[0],self.N[1],self.N[2]):
                ivec = self.to_BZ(np.array([i,j,k]))
                if not self.is_wave_cancel(ivec,g):
                    n += 1
            Nw = n
            w = np.zeros((self.dim,Nw))
            G2 = np.zeros(Nw)
            n = 0
            for (i,j,k) in np.ndindex(self.N[0],self.N[1],self.N[2]):
                ivec = self.to_BZ(np.array([i,j,k]))
                if not self.is_wave_cancel(ivec,g):
                    w[:,n] = ivec
                    G2[n] = wave_norm(ivec,self.shape)
                    n += 1

        # Sort G2 in ascending order, returned the corresponding indices
        ind = np.argsort(G2)
        self.Nw = Nw
        self.waves = w[:,ind]
        self.Gsq = G2[ind]

def wave_norm(G,shape):
    """
    G is a wave vector with 1, 2, or 3 elements
    shape is a shape matrix
    """
    v = np.dot(G,shape.g)
    return np.dot(v,v)

