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

2. Usage
^^^^^^^^

::

    >>>import gyroid as gy

    >>>uc = gy.UnitCell(3)
    >>>group = gy.Group(3,gy.BRAVAIS,uc.shape,"Ia-3d")
    >>>grid = gy.Grid(np.array([4,4,4]),group)
    >>>basis = gy.Basis(group,grid)

    >>>render_structure_3d(basis,32,32,32,1.0)

Ask for Help
------------

* You can directly contact me at liuyxpp@gmail.com.
* You can join the mailinglist by sending an email to gyroid@librelist.com and replying to the confirmation mail. To unsubscribe, send a mail to gyroid-unsubscribe@librelist.com and reply to the confirmation mail.

Links
-----

* `Website <http://liuyxpp.bitbucket.org>`_
* `Development version <http://bitbucket.org/liuyxpp/gyroid/>`_

'''
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='gyroid',
    version='0.2',
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

