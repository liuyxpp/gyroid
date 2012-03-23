# -*- coding: utf-8 -*-
"""
gyroid.util
===============

"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import colors
from mayavi import mlab

__all__ = ["render_structure_1d","render_structure_2d","render_structure_3d"]

def render_structure_1d(basis,Na,c,
                        save_data=True,data_name="struct1d.mat",
                        save_img=True,show_img=True,
                        img_name="struct1d.png",
                        **kwargs):
    """
    basis - An instance of Basis class
    Na    - Number of grids to discretes the side a of the unit cell.
    c     - coefficients for each basis function.
    """

    struct = basis.generate_structure(Na,c)
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

def render_structure_2d(basis,Na,Nb,c,
                        save_data=True,data_name="struct2d.mat",
                        save_img=True,show_img=True,
                        img_name="struct2d.png",
                        levels=None,cmap=None,
                        **kwargs):
    """
    basis - An instance of Basis class
    Na,Nb - Number of grids to discrete the side a of the unit cell.
    c     - coefficients for each basis function.
    """

    struct = basis.generate_structure((Na,Nb),c)
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
    if save_img:
        plt.savefig(img_name)
    if show_img:
        plt.show()
    return rx,ry,struct

def render_structure_3d(basis,Na,Nb,Nc,c,
                        save_data=True,data_name="struct3d.mat",
                        save_img=True,show_img=True,
                        img_name="struct3d.png",
                        levels=None,cmap=None,
                        **kwargs):
    """
    basis    - An instance of Basis class
    Na,Nb,Nc - Number of grids to discrete the side a of the unit cell.
    c        - coefficients for each basis function.
    NOTE: the best way to view 3D volume data is: first save the data to mat, and let Matlab (C) render the volume data.
    """

    struct = basis.generate_structure((Na,Nb,Nc),c)
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

