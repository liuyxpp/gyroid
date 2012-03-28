# -*- coding: utf-8 -*-
"""
gyroid.group
============

"""

import numpy as np
from numpy.linalg import inv

#from .symmetry import find_symmetry_for_group
from .common import HM_ITA_TABLE1,HM_ITA_TABLE2,HM_ITA_TABLE3
from .common import GROUP_ORDER_TABLE1,GROUP_ORDER_TABLE2,GROUP_ORDER_TABLE3
from .space_group import symmetry_generator1,symmetry_generator2,symmetry_generator3

__all__ = ["Group","CayleyTable"]

class Group(object):
    ''' A :class:`Group` object is a representaion of a space group.

    All symmetries in a space group must have the same basis type, i.e. they
    must all be either the BRAVAIS or the CARTESIAN bases.

    The space group is constructed either by providing a Hermann-Mauguin Symbol (HM_symbol) or a sequential number as given in the International Tables for Crystallography, Vol. A (ITA_number)

    :members: dim, type, shape, symm

    '''

    def __init__(self,dim,basis_type,shape_matrix,
                 HM_symbol=None,ITA_number=None):
        '''
        If both HM symbol and ITA number are provided, use HM symbol only. And the consistency of HM symbol and ITA number is checked.
        '''

        if dim < 1 or dim > 3:
            raise ValueError('Only 1D, 2D, and 3D space allowed.')
        self.dim = dim
        self.type = basis_type
        self.shape = shape_matrix
        if HM_symbol:
            self.__create_group_by_symbol(HM_symbol)
            return
        if ITA_number is None:
            return
        self.__create_group_by_number(ITA_number)

    def __find_ITA_number(self,HM_symbol):
        if self.dim == 1 and HM_symbol in HM_ITA_TABLE1:
            return HM_ITA_TABLE1.index(HM_symbol) + 1
        if self.dim == 2 and HM_symbol in HM_ITA_TABLE2:
            return HM_ITA_TABLE2.index(HM_symbol) + 1
        if self.dim == 3 and HM_symbol in HM_ITA_TABLE3:
            return HM_ITA_TABLE3.index(HM_symbol) + 1
        return 0

    def __is_valid_ITA_number(self,ITA_number):
        if self.dim == 1 and ITA_number > 0 and ITA_number < 3:
            return True
        if self.dim == 2 and ITA_number > 0 and ITA_number < 18:
            return True
        if self.dim == 3 and ITA_number > 0 and ITA_number < 231:
            return True
        return False

    def __find_HM_symbol(self,ITA_number):
        if self.dim == 1:
            return HM_ITA_TABLE1[ITA_number-1]
        if self.dim == 2:
            return HM_ITA_TABLE2[ITA_number-1]
        if self.dim == 3:
            return HM_ITA_TABLE3[ITA_number-1]

    def __find_group_order(self,ITA_number):
        if self.dim == 1:
            return GROUP_ORDER_TABLE1[ITA_number-1]
        if self.dim == 2:
            return GROUP_ORDER_TABLE2[ITA_number-1]
        if self.dim == 3:
            return GROUP_ORDER_TABLE3[ITA_number-1]

    def __find_symmetry_for_group(self,ITA_number):
        if self.dim == 1:
            return symmetry_generator1(ITA_number,self.type,self.shape)
        if self.dim == 2:
            symm = symmetry_generator2(ITA_number,self.type,self.shape)
        if self.dim == 3:
            symm = symmetry_generator3(ITA_number,self.type,self.shape)

        # generating group symmetries using generators
        cayley = CayleyTable(symm)
        for iterate in range(20):
            is_group = True
            # Add missing inverses of existing elements to group
            k = len(symm)
            for i in range(k):
                isymm = symm[i].inverse()
                if isymm not in symm:
                    symm.append(isymm)
                    is_group = False
            # Add products of existing elements to group
            cayley.update(symm)
            k = len(symm)
            for (i,j) in np.ndindex(k,k):
                #if (cayley.index[i,j] == -1) and (cayley.table[i][j] not in symm):
                if cayley.table[i][j] not in symm:
                    symm.append(cayley.table[i][j])
                    is_group = False
            # No more symmetry element is added before end of iteration.
            if is_group:
                return symm
        raise ValueError("Failed to create a group.")

    def __create_group_by_symbol(self,HM_symbol):
        number = self.__find_ITA_number(HM_symbol)
        self.__create_group_by_number(number)

    def __create_group_by_number(self,ITA_number):
        if not self.__is_valid_ITA_number(ITA_number):
            raise ValueError("Incorrect HM symbol or ITA number.")
        self.number = ITA_number
        self.HM_symbol = self.__find_HM_symbol(ITA_number)
        self.symm = self.__find_symmetry_for_group(ITA_number)
        self.order = len(self.symm)
        if self.order != self.__find_group_order(ITA_number):
            raise ValueError("Incorrect order for group")


class CayleyTable(object):
    """ Make a full Cayley table from a list of :class:`Symmetry` objects.

    """

    def __init__(self,symm=None):
        # create a self.order x self.order table
        # NOTE: table is not an numpy ndarray!
        if symm is None:
            self.table = [[]]
        else:
            self.order = len(symm)
            self.table = [[symm[i]*symm[j] for i in range(self.order)]
                                           for j in range(self.order)]

    def update(self,symm):
        ''' Update the :class:`CayleyTable` object.

        Currently, `len(symm)` must be larger than `self.order`.

        '''

        l0 = self.order
        l1 = len(symm)
        if l1 <= l0:
            return

        # Extend original rows
        for i in range(l0):
            self.table[i].extend([symm[i]*symm[j] for j in range(l0,l1)])
        # Add additional rows
        self.table.extend([[symm[i]*symm[j] for j in range(l1)]
                                            for i in range(l0,l1)])
        self.order = l1


