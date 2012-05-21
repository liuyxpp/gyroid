'''
gyroid
======

**gyroid** is a python package that generates *symmetry adapted basis functions* based on the space group of a unit cell.

Quickstart
----------

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

Ask for Help
------------

* You can directly contact me at liuyxpp@gmail.com.
* You can join the mailinglist by sending an email to gyroid@librelist.com and replying to the confirmation mail. To unsubscribe, send a mail to gyroid-unsubscribe@librelist.com and reply to the confirmation mail.

Links
-----

* `Documentation <http://packages.python.org/gyroid>`_
* `Website <http://liuyxpp.bitbucket.org>`_
* `Development version <http://bitbucket.org/liuyxpp/gyroid/>`_

'''
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='gyroid',
    version='0.4',
    license='BSD',
    description='A symmetry adapted basis function (SABF) generator.',
    author='Yi-Xin Liu',
    author_email='liuyxpp@gmail.com',
    url='https://bitbucket.org/liuyxpp/gyroid',
    packages=['gyroid'],
    include_package_data=True,
    zip_safe=False,
    long_description=__doc__,
    platform='linux',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib>=1.0.1',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Education',
        'Topic :: Multimedia :: Graphics',
    ]
     )

