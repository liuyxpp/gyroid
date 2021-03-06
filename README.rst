gyroid
======

**gyroid** is a python package that generates *symmetry adapted basis functions* based on the space group of a unit cell.

Quickstart
----------

1. Install
^^^^^^^^^^

::

    $ pip install gyroid

or

::

    $ tar -xvf gyroid-xxx.tar.gz
    $ cd gyroid-xxx
    $ python setup.py install

Required packages:

* `numpy`: it should be installed before installing gyroid.
* `scipy`: use it to save data in Matlab mat format.
* `matplotlib`: 2D Graphic plotting.
* `mayavi`: it depends on many packages, e.g. VTK (compiled with python wrapper and shared library on). If you do not need the `render_structure_3d` functionality, simply ignore it.

2. Usage
^^^^^^^^

Following is a typical usange of the package. It will produce a set of SABFs with point group Ia-3d (#230) in a cubic unit cell. The last line will syntheses a gyroid structure with all coefficients for SABF equal to 1.0, save the structure data into a Matlab mat file, show a screenshot of the rendered image and save the image in a file.

>>> import gyroid as gy
>>> import numpy as np
>>> N1,N2,N3 = 32,32,32
number of grids in each dimension of a unit cell
>>> uc = gy.UnitCell(3)
create a standard cubic unit cell with side length 1.0
>>> group = gy.Group(3,gy.BRAVAIS,uc.shape,"Ia-3d")
create a Ia-3d point group
>>> grid = gy.Grid(np.array([N1,N2,N3]),group)
create a collection of waves that are not canceled
>>> basis = gy.Basis(group,grid)
create the SABFs
>>> gy.render_structure_3d(basis,grid,N1,N2,N3,1.0)
create the gyroid structure and render it

Ask for Help
------------

* You can directly contact me at lyx@fudan.edu.cn.
* You can file an issue at github.com.

Links
-----

* `Documentation <http://packages.python.org/gyroid>`_
* `Development version <https://github.com/liuyxpp/gyroid>`_

