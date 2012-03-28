
Welcome to gyroid
=================

Welcome to gyroid's documentation.

About gyroid
------------

`gyroid` is a python package that generates *symmetry adapted basis functions* (SABF) based on the point group of a unit cell.

It was originally implemented in the spirit of creating ordered structures from SABF as the initial guess for polymer self-consistent field theory (SCFT) calculations. Many of the algorithms are adopted from the software package `pscf <http://www.cems.umn.edu/research/morse/code/pscf/home.php>`_ authored by David Morse, etc. written in Fortran.

Quick Start
-----------

1. Install
^^^^^^^^^^

::

    $ easy_install gyroid

or

::

    $ tar -xvf gyroid-xxx.tar.gz
    $ cd gyroid-xxx
    $ python setup.py install

Required packages:

* `numpy`: it should be installed before installing gyroid.
* `scipy`: use it to save data in Matlab mat format.
* `matplotlib`: 2D Graphic plotting.
* `mayavi`: it depends on many packages, e.g. VTK (compiled with python wrapper and shared library on). If you do not need the render_structure_3d functionality, simply ignore it.

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

API Reference
-------------

For advanced usage of this package, or if you want to modify the package yourself (e.g. you can add unsupported piont group into the package by modifying `gyroid.space_group.py` module.), please refer to the full API documentation.

.. toctree::
   :maxdepth: 2

   gyroid

Indices
-------

* :ref:`genindex`
* :ref:`search`

Ask for Help
------------

* You can directly contact me at liuyxpp@gmail.com.
* You can join the mailinglist by sending an email to gyroid@librelist.com and replying to the confirmation mail. To unsubscribe, send a mail to gyroid-unsubscribe@librelist.com and reply to the confirmation mail.

Links
-----

* `Documentation <http://packages.python.org/gyroid>`_
* `Website <http://liuyxpp.bitbucket.org>`_
* `Development version <http://bitbucket.org/liuyxpp/gyroid/>`_

