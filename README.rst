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

