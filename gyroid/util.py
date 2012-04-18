# -*- coding: utf-8 -*-
"""
gyroid.util
===========

"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import colors
from mayavi import mlab

from .unitcell import UnitCell
from .group import Group
from .grid import Grid
from .basis import Basis

__all__ = [
    "render_structure_1d",
    "render_structure_2d",
    "render_structure_3d",
    "prepare_scft_input"]

def prepare_scft_input(dim,grid_num_vec,cryst_system,
                       cryst_param_vec,sym_group,basis_grid_vec,basis_c,
                       data_file="field_in.mat",show_img=False,
                       save_img=False,img_file="field_in.png",
                       **kwargs):
    b = "Bravais"
    uc = UnitCell(dim,cryst_system,cryst_param_vec);
    g = Group(dim,b,uc.shape,sym_group)
    gd = Grid(basis_grid_vec,g)
    bs = Basis(g,gd)

    c = np.zeros(bs.N)
    N = basis_c.size
    if N < bs.N:
        c[0:N] = basis_c
    else:
        c = basis_c[0:bs.N]

    if dim == 1:
        render_structure_1d(bs,gd,grid_num_vec[0],c,
                            data_name=data_file,save_img=save_img,
                            show_img=show_img,img_name=img_file,
                            **kwargs)
    if dim == 2:
        render_structure_2d(bs,gd,grid_num_vec[0],grid_num_vec[1],c,
                            data_name=data_file,save_img=save_img,
                            show_img=show_img,img_name=img_file,
                            **kwargs)
    if dim == 3:
        render_structure_3d(bs,gd,grid_num_vec[0],grid_num_vec[1],
                            grid_num_vec[2],c,
                            data_name=data_file,save_img=save_img,
                            show_img=show_img,img_name=img_file,
                            **kwargs)

def render_structure_1d(basis,grid,Na,c,
                        save_data=True,data_name="struct1d.mat",
                        save_img=True,show_img=True,
                        img_name="struct1d.png",
                        **kwargs):
    ''' Calculate and render 1D structure for given SABF and unit cell.

    :param basis: a set of SABFs
    :type basis: :class:`Basis`
    :param Na: number of grids in **a** of the unit cell.
    :type Na: integer
    :param c: coefficients for each SABF
    :type c: 1D `numpy.array`
    :param save_data: if True, save data in file with Matlab mat format
    :type save_data: bool
    :param data_name: the file name of the data file
    :type data_name: string
    :param save_img: if True, save image in file, the format is determined by the extension of the image file name
    :type save_img: bool
    :param img_name: the file name of the image file
    :type img_name: string
    :param show_img: if True, show image on the screen
    :type show_img: bool
    :param kwargs: any extra key words arguments will be passed to plot functions

    '''

    #struct = basis.generate_structure(Na,c)
    struct = basis.generate_structure_by_fft((Na,),c,grid)
    # For debug only
    #print basis.fft2sabf(np.fft.fftn(struct),grid)
    a = 1.0 * basis.shape.h[0,0]
    rx = np.array([a*i/Na for i in np.arange(Na)])

    if save_data:
        scipy.io.savemat(data_name,{"rx":rx,"struct":struct})
    if save_img or show_img:
        plt.plot(rx,struct,**kwargs)
    if save_img:
        plt.savefig(img_name)
    if show_img:
        plt.show()
    return rx,struct

def render_structure_2d(basis,grid,Na,Nb,c,
                        save_data=True,data_name="struct2d.mat",
                        save_img=True,show_img=True,
                        img_name="struct2d.png",
                        levels=None,cmap=None,
                        **kwargs):
    ''' Calculate and render 2D structure for given SABF and unit cell.

    :param basis: a set of SABFs
    :type basis: :class:`Basis`
    :param Na: number of grids in **a** of the unit cell.
    :type Na: integer
    :param Nb: number of grids in **b** of the unit cell.
    :type Nb: integer
    :param c: coefficients for each SABF
    :type c: 1D `numpy.array`
    :param save_data: if True, save data in file with Matlab mat format
    :type save_data: bool
    :param data_name: the file name of the data file
    :type data_name: string
    :param save_img: if True, save image in file, the format is determined by the extension of the image file name
    :type save_img: bool
    :param img_name: the file name of the image file
    :type img_name: string
    :param show_img: if True, show image on the screen
    :type show_img: bool
    :param kwargs: any extra key words arguments will be passed to plot functions

    '''

    # If generate_structure_by_fft failed
    # Give generate_structure a try.
    #struct = basis.generate_structure((Na,Nb),c)
    struct = basis.generate_structure_by_fft((Na,Nb),c,grid)

    # For debug only
    print "Input c: ",c
    print "c from constructed structure: "
    print basis.fft2sabf(np.fft.fftn(struct),grid)

    rx = np.zeros((Na,Nb))
    ry = np.zeros((Na,Nb))
    for (i,j) in np.ndindex(Na,Nb):
        x = (1.0*np.array([i,j])) / (Na,Nb)
        rx[i,j],ry[i,j] = np.dot(x,basis.shape.h)

    if save_data:
        scipy.io.savemat(data_name,{"rx":rx,"ry":ry,"struct":struct})
    if save_img or show_img:
        dx = rx.max() - rx.min()
        dy = ry.max() - ry.min()
        w,h = plt.figaspect(float(dy/dx)) # float is must
        # No frame, white background, w/h aspect ratio figure
        fig = plt.figure(figsize=(w,h),frameon=False,
                         dpi=80,facecolor='w')
        # full figure subplot, no border, no axes
        ax = fig.add_axes([0,0,1,1],frameon=False,axisbg='w')
        # no ticks
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        # Default: there are 256 contour levels
        if levels is None:
            step = (struct.max() - struct.min()) / 256
            levels = np.arange(struct.min(),struct.max()+step,step)
        # Default: colormap is monochromatic red
        if cmap is None:
            clr = np.zeros((256,3))
            for i in np.arange(256):
                clr[i,0] = i / 255.0
            cmap = colors.ListedColormap(clr)
        # actual plot
        ax.contourf(rx,ry,struct,levels=levels,
                    cmap=cmap,antialiased=False,**kwargs)
        #ax.contourf(rx,ry,struct)
    if save_img:
        plt.savefig(img_name)
    if show_img:
        plt.show()
    return rx,ry,struct

def render_structure_3d(basis,grid,Na,Nb,Nc,c,
                        save_data=True,data_name="struct3d.mat",
                        save_img=True,show_img=True,
                        img_name="struct3d.png",
                        levels=None,cmap=None,
                        **kwargs):
    ''' Calculate and render 3D structure for given SABF and unit cell.

    :param basis: a set of SABFs
    :type basis: :class:`Basis`
    :param Na: number of grids in **a** of the unit cell.
    :type Na: integer
    :param Nb: number of grids in **b** of the unit cell.
    :type Nb: integer
    :param Nc: number of grids in **c** of the unit cell.
    :type Nc: integer
    :param c: coefficients for each SABF
    :type c: 1D `numpy.array`
    :param save_data: if True, save data in file with Matlab mat format
    :type save_data: bool
    :param data_name: the file name of the data file
    :type data_name: string
    :param save_img: if True, save image in file, the format is determined by the extension of the image file name
    :type save_img: bool
    :param img_name: the file name of the image file
    :type img_name: string
    :param show_img: if True, show image on the screen
    :type show_img: bool
    :param kwargs: any extra key words arguments will be passed to plot functions

    '''

    #struct = basis.generate_structure((Na,Nb,Nc),c)
    struct = basis.generate_structure_by_fft((Na,Nb,Nc),c,grid)
    # For debug only
    #print basis.fft2sabf(np.fft.fftn(struct),grid)
    rx = np.zeros((Na,Nb,Nc))
    ry = np.zeros((Na,Nb,Nc))
    rz = np.zeros((Na,Nb,Nc))
    for (i,j,k) in np.ndindex(Na,Nb,Nc):
        x = (1.0*np.array([i,j,k])) / (Na,Nb,Nc)
        rx[i,j,k],ry[i,j,k],rz[i,j,k] = np.dot(x,basis.shape.h)

    if save_data:
        scipy.io.savemat(data_name,
                         {"rx":rx,"ry":ry,"rz":rz,"struct":struct})
    if save_img or show_img:
        mlab.contour3d(rx,ry,rz,struct,**kwargs)
    if save_img:
        mlab.savefig(img_name)
    if show_img:
        plt.show()
    return rx,ry,rz,struct

