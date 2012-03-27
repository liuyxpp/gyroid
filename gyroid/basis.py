# -*- coding: utf-8 -*-
"""
gyroid.basis
===============

:copyright: (c) 2012 by Yi-Xin Liu
:license: BSD, see LICENSE for more details.

"""

import numpy as np
from numpy.linalg import inv

from .common import BRAVAIS,CARTESIAN,EPS

__all__ = ["Basis","StarSet","StarAtom"]

class Basis(object):
    """
    """
    def __init__(self,group,grid):
        self.dim = group.dim
        self.shape = group.shape
        self.stars = []
        self.starmap = {} # key: G (within BZ), value: index of a star
        self.N = 0  # this is for coefficients, N = #(closed star) +
                    # 2 * #(open star pair)
        n = 0   # record the number of StarAtom
        m = 0   # record the index of coefficent
        G2_pre = grid.Gsq[0] # Previous G2
        for G2 in grid.Gsq:
            if np.abs(G2-G2_pre) > EPS:
                s = StarSet(group,grid,G2_pre)
                # s.stars is an Python list.
                self.stars.extend(s.stars)
                # map G to the index of star
                i = 0
                for star in s.stars:
                    iw = 0
                    for G in star.waves.T:
                        key = tuple(G.astype(int))
                        self.starmap[key] = (n+i, iw, 0)
                        iw += 1
                    if star.iwaves is not None:
                        iw = 0
                        for Gi in star.iwaves.T:
                            key = tuple(Gi.astype(int))
                            self.starmap[key] = (n+i, iw, 1)
                            iw += 1
                    i += 1
                n += i  # i stars has been counted in last cycle
                self.N += s.N
                G2_pre = G2
#        i,iw,flag = index_stars((0,0,-15),self.stars)
#        print i,iw,flag
#        print self.stars[i].waves
#        print self.stars[i].iwaves

    def generate_structure(self,real_grid,c):
        """
        real_grid is a tuple contains the discrete numbers for a, b, and c unit vectors in real space.
        c is an array contains the coefficients for each basis.

        NOTE: this method is much slower than generate_structure_by_fft()
        """
        if np.size(real_grid) != self.dim:
            raise ValueError("Dimension not match when generating struct.")
        if np.size(c) == 1:
           cc = c
           c = np.zeros(self.N)
           c.fill(cc)
        elif np.size(c) != self.N:
            raise ValueError("No. of Bases and coefficients not match when generating structure.")

        struct = np.zeros(real_grid)
        for ind,v in np.ndenumerate(struct):
            # numpy.ndarray can divide tuple
            #x = ind * self.shape.l / real_grid
            x = 1.0 * np.array(ind) / real_grid
            i = 0
            for s in self.stars:
                if s.iwaves is None:
                    f1,f2 = s.f(x,self.shape,c[i],c[i])
                else:
                    f1,f2 = s.f(x,self.shape,c[i],c[i+1])
                struct[ind] += f1
                i += 1
                if f2 is not None:
                    struct[ind] += f2
                    i += 1
        vol = np.dot(struct.shape,struct.shape)
        return struct/vol

    def generate_structure_by_fft(self,real_grid,c,grid):
        if np.size(real_grid) != self.dim:
            raise ValueError("Dimension not match when generating struct.")
        if np.size(c) == 1:
           cc = c
           c = np.zeros(self.N)
           c.fill(cc)
        elif np.size(c) != self.N:
            raise ValueError("No. of Bases and coefficients not match when generating structure.")

        c_fft = self.sabf2fft(c,real_grid,grid)
        return np.fft.ifftn(c_fft).real

    def sabf2fft(self,c,fft_grid,grid):
        """
        c           - is a list of coefficients for SABF.
        fft_grid    - tuple with discretized number for each side.
        """

        if np.size(c) != self.N:
            raise ValueError("No. of input coefficients != No. of stars.")
        if np.size(fft_grid) != self.dim:
            raise ValueError("Dimension not match in sabf2fft.")

        sqr2 = np.sqrt(2.0)
        c_fft = np.zeros(fft_grid).astype(complex)
        for ind in np.ndindex(fft_grid):
            G = np.array(ind)
            G,G2 = grid.to_BZ(G)
            key = tuple(G.astype(int))
            if self.starmap.has_key(key):
                i, iw, flag = self.starmap[key] #index_stars(G,self.stars)
            else:
                i, iw, flag = None,None,None
            if i is not None:
                if flag == 0:
                    if self.stars[i].iwaves is None:
                        c_fft[ind] = self.stars[i].c[iw] * c[i]
                    else:
                        c_fft[ind] = self.stars[i].c[iw] * (complex(
                                                c[i],-c[i+1]) / sqr2)
                        print i,iw,flag
                else:
                    c_fft[ind] = self.stars[i].ic[iw] * (complex(
                                                c[i-1],c[i]) / sqr2)
        return c_fft

    def fft2sabf(self,c_fft,grid):
        if np.ndim(c_fft) != self.dim:
            raise ValueError("Dimension not match in sabf2fft.")

        fft_grid = np.shape(c_fft)
        sqr2 = np.sqrt(2.0)
        c = np.zeros(self.N)
        i = 0
        for s in self.stars:
            G = s.waves.T[0]
            if G[0] > fft_grid[0]/2:
                G = (-G) % fft_grid
                ind = tuple(G.astype(int))
                z = c_fft[ind].conjugate()
            else:
                ind = tuple(G.astype(int))
                z = c_fft[ind]

            if s.iwaves is None:
                c[i] = (z/s.c[0]).real
                i += 1
            else:
                c[i] = sqr2 * (z/s.c[0]).real
                Gi = s.iwaves.T[s.N-1]
                if Gi[0] > fft_grid[0]/2:
                    Gi = (-Gi) % fft_grid
                    ind = tuple(Gi.astype(int))
                    z = c_fft[ind].conjugate()
                else:
                    ind = tuple(Gi.astype(int))
                    z = c_fft[ind]
                c[i+1] = sqr2 * (z/s.ic[s.N-1]).imag
                i += 2
        return c


class StarSet(object):
    """
    A StarSet is a collection of stars with all waves have same magnitude.
    For example, for P6mm in a 32 x 32 grid HEXAGONAL unit cell, an example of a
    collection of star with same magnitude is:
        [[ 8  8 7  7 5  5 3  3 0  0 -3 -3 -5 -5 -7 -7 -8 -8],
         [-3 -5 0 -7 3 -8 5 -8 7 -7  8 -5  8 -3  7  0  5  3]]
    It has two *closed* stars:
        [[ 8  8 5  5 3  3 -3 -3 -5 -5 -8 -8],
         [-3 -5 3 -8 5 -8  8 -5  8 -3  5  3]]
    and
        [[7  7 0  0 -7 -7],
         [0 -7 7 -7  7  0]]
    """

    def __init__(self,group,grid,Gsq):
        if self.__check_cancel():
            raise ValueError("Check cancel failed when creating a Star.")
        self.dim = group.dim
        self.Gsq = Gsq
        waves = self.__select_waves(grid,Gsq)
        sorted_waves,phases = self.__sort_waves(waves)
        self.__find_stars(group,grid,sorted_waves)

    def __select_waves(self,grid,G2):
        (ind,) = np.where(np.abs(grid.Gsq-G2)<EPS)
        if np.max(ind) - np.min(ind) + 1 != np.size(ind):
            raise ValueError("Waves in Grid not sorted according to G^2.")
        return grid.waves[:,ind]

    def __check_cancel(self):
        """
        Not implemented yet. Since we have excluded the cancel waves in
        creating Grid.waves.
        Currently, we do not support canceled stars.
        """
        return False

    def __sort_waves(self,waves,phases=None):
        if self.dim == 1:
            if phases is None:
                if np.size(waves,1) == 1:
                    return waves,None
                return (np.fliplr(np.sort(waves)),None)
            else:
                pw = np.vstack([phases,waves])
                ind = np.lexsort(pw)
                pw_sorted = np.fliplr(pw.take(ind,axis=-1))
                return (np.array([pw_sorted[1]]),pw_sorted[0])

        if self.dim == 2:
            if phases is None:
                rw = np.vstack([waves[1],waves[0]])
                ind = np.lexsort(rw)
                return (np.fliplr(waves.take(ind,axis=-1)),None)
            else:
                prw = np.vstack([phases,waves[1],waves[0]])
                ind = np.lexsort(prw)
                prw_sorted = np.fliplr(prw.take(ind,axis=-1))
                return (np.vstack([prw_sorted[2],prw_sorted[1]]),
                        prw_sorted[0])

        if self.dim == 3:
            if phases is None:
                rw = np.vstack([waves[2],waves[1],waves[0]])
                ind = np.lexsort(rw)
                return (np.fliplr(waves.take(ind,axis=-1)),None)
            else:
                prw = np.vstack([phases,waves[2],waves[1],waves[0]])
                ind = np.lexsort(prw)
                prw_sorted = np.fliplr(prw.take(ind,axis=-1))
                return (np.vstack([
                    prw_sorted[3],prw_sorted[2],prw_sorted[1]]),
                    prw_sorted[0])

        # Following code is a shortcut but hard to read
        # ind = np.lexsort(waves.T)
        # return np.fliplr(np.fliplr(waves.T.take(ind,axis=-1)).T)

    def __calc_phase(self,G,t,basis_type):
        twopi = 2.0 * np.pi
        if basis_type == BRAVAIS:
            return twopi * np.round(np.dot(G,t)).astype(type(G[0]))
        else:
            return np.dot(G,t)

    def __calc_wave(self,G,R,basis_type):
        if basis_type == BRAVAIS:
            return np.round(np.dot(G,R)).astype(type(G[0]))
        else:
            return np.dot(G,R)

    def __form_star(self,G,group,grid,waves):
        star_waves = None
        phases = None
        for i in np.arange(group.order):
            Gn = self.__calc_wave(G,group.symm[i].R,group.type)
            # Pseudo-Spectral method
            Gn,Gn2 = grid.to_BZ(Gn)

            if index_waves(Gn,waves.T) is not None:
                if star_waves is None:
                    star_waves = np.array([Gn])
                    ph = self.__calc_phase(G,group.symm[i].t,group.type)
                    phases = np.array([ph])
                else:
                    if index_waves(Gn,star_waves) is None:
                        star_waves = np.append(star_waves,[Gn],axis=0)
                        ph = self.__calc_phase(G,group.symm[i].t,group.type)
                        phases = np.append(phases,[ph],axis=0)
            else:
                raise ValueError("Waves does not contain entire star.")
        return star_waves.T,phases

    def __find_stars(self,g,grid,waves):
        """
        For waves with a same |G|^2, they may form a closed star, two open
        stars, or several closed stars.
        """

        self.stars = []
        self.N = 0
        rw = waves
        while rw is not None:
            G1 = rw[:,0]
            star_waves, phases = self.__form_star(G1,g,grid,rw)
            star_waves, phases = self.__sort_waves(star_waves,phases)
            Gi = -1.0 * G1
            Gi,Gi2 = grid.to_BZ(Gi)
            if index_waves(Gi,star_waves.T) is not None:
                # a closed star
                self.stars.append(StarAtom(grid,self.Gsq,star_waves,phases))
                self.N += 1
                if np.size(rw,1) == np.size(star_waves,1):
                    return
                tw = None
                for w in rw.T:
                    if index_waves(w,star_waves.T) is None:
                        if tw is None:
                            tw = np.array([w])
                        else:
                            tw = np.append(tw,[w],axis=0)
                rw = tw.T
            else:
                # an open star pair
                invert_waves, invert_phases = self.__form_star(
                                                    Gi,g,grid,rw)
                invert_waves, invert_phases = self.__sort_waves(
                                            invert_waves,invert_phases)
                self.stars.append(StarAtom(grid,self.Gsq,
                    star_waves,phases,invert_waves,invert_phases))
                self.N += 2
                if np.size(rw,1) == np.size(star_waves,1) + np.size(
                                                        invert_waves,1):
                    return
                tw = None
                for w in rw.T:
                    if index_waves(w,star_waves.T) is None and index_waves(w,invert_waves.T) is None:
                        if tw is None:
                            tw = np.array([w])
                        else:
                            tw = np.append(tw,[w],axis=0)
                rw = tw.T


class StarAtom(object):
    """
    A StarAtom is a star or an open star pair.
    """

    def __init__(self,grid,Gsq,waves,phases,iwaves=None,iphases=None):
        self.N = np.size(waves,1)
        self.Gsq = Gsq
        if iwaves is not None:
            if iphases is None:
                raise ValueError(
                    "Coefficients expected for inverted star.")
            if np.size(iwaves,1) != self.N:
                raise ValueError(
                    "No. Waves in star and inverted star not match.")
        else:
            if iphases is not None:
                raise ValueError(
                    "Waves expected for inverted coefficients.")
        self.waves = waves
        self.iwaves = iwaves
        self.c, self.ic = self.__find_coeff(phases,iphases)
        if self.iwaves is None:
            self.__set_coeff_for_closed_star(grid)

    def f(self,x,shape,c1,c2):
        """
        the value of basis function f(r) at position r.
        x   - a Bravais type real space vector.
        c1  - the coefficient for waves
        c2  - the coefficient for iwaves
        """
        if self.iwaves is None:
            f1 = self.__f(self.waves,self.c,x,shape)
            return (c1*f1.real,None)
        else:
            v1 = self.__f(self.waves,self.c,x,shape)
            v2 = self.__f(self.iwaves,self.ic,x,shape)
            f1 = (v1 + v2) / np.sqrt(2.0)
            f2 = complex(0.0,1.0) * (v1 - v2) / np.sqrt(2.0)
            return (c1*f1.real,c2*f2.real)

    def __f(self,waves,c,x,shape):
        f = 0
        i = 0
        for G in waves.T:
            gr = np.dot(np.dot(G,shape.g),np.dot(x,shape.h))
            f += c[i] * np.exp(complex(0.0,1.0) * gr)
            i += 1
        return f

    def __find_coeff(self,phases,iphases):
        if iphases is None:
            # c_norm = exp(i*phi)
            c_norm = np.exp(complex(0.0,phases[0])) * np.sqrt(self.N)
            return ([np.exp(complex(0.0,phases[i]))/c_norm
                    for i in np.arange(self.N)], None)
        else:
            c_norm = np.exp(complex(0.0,phases[0])) * np.sqrt(self.N)
            ic_norm = np.exp(complex(0.0,phases[self.N-1]))*np.sqrt(self.N)
            return ([np.exp(complex(0.0,phases[i]))/c_norm
                    for i in np.arange(self.N)],
                    [np.exp(complex(0.0,iphases[i]))/ic_norm
                    for i in np.arange(self.N)]
                   )

    def __set_coeff_for_closed_star(self,grid):
        """
        For an ordered closed star, if we denote the first wave in the star G1, then its inversion (-G1) must be the last wave in the star.
        """
        G = self.waves[:,0]
        Gi = -1.0 * G
        Gi,Gi2 = grid.to_BZ(Gi)
        # find index of Gi in the star
        i = index_waves(Gi,self.waves.T)

        c, ci = self.c[0], self.c[i]
        if np.abs(c.imag) < EPS:
            c1 = c.real
        else:
            raise ValueError("""First coefficient in closed star has
                              imaginary part.""")
        if np.abs(ci.imag) < EPS:
            c2 = ci.real
        else:
            raise ValueError("""Last coefficient in closed star has
                              imaginary part.""")

        # for inversion star pairs, the first star's sign is +1,
        # the next is -1.
        if np.abs(c1 - c2) < EPS:
            return 1
        elif np.abs(c1 + c2) < EPS:
            self.c = self.c * complex(0.0,-1.0)
            return -1
        else:
            raise ValueError("""Closed star is neither cosine-like nor
                             sine-like.""")

def index_waves(w,waves):
    """
    w is a row vector.
    waves is a 2D array, each row is a row vector.
    """
    if np.size(w) != np.size(waves,1):
        return None
    if waves is None:
        return None
    i = 0
    for ww in waves:
        if np.all(np.abs(ww-w) < EPS):
            return i
        i += 1
    return None

def index_stars(G,stars):
    """
    """
    if np.size(G) != np.size(stars[0].waves.T[0]):
        return None,None,None
    i = 0
    for s in stars:
        iw = index_waves(G,s.waves.T)
        if iw is not None:
            return i,iw,0
        if s.iwaves is not None:
            iw = index_waves(G,s.iwaves.T)
            if iw is not None:
                return i,iw,1
        i += 1
    return None,None,None

