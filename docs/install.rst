.. _install:

Installation
============

lmfit2 is availble in both C and in Python, with each having a different
installation procedure. The C version needs to be compiled against the
`RST package <https://github.com/SuperDARN/rst>`_. The Python version is
not a package, but instead a python script.

The C version is the recommended version to be used. The Python version
is only provided for testing purposes.


C Version Installation
----------------------

First, `install RST using this guide <https://radar-software-toolkit-rst.readthedocs.io/en/latest/>`_.

Next, obtain a copy of lmfit2 source code. Here we'll assume you are
using `git`, then compile the C source code:

.. code-block:: bash

    BASEPATH=`pwd`
    git clone https://github.com/asreimer/lmfit2.git@v1.0
    cd ${BASEPATH}/lmfit2/C
    LIBPATH=$RSTPATH/lib IPATH=$RSTPATH/include make

After the build completes, a `make_lmfit2` binary will be available in
`${BASEPATH}/lmfit2/C/bin`. You can add this location to your user path or you
can copy the binary to an existing location on your path.


Test the C binary
-----------------

To make sure that the installation went alright, you can execute a unit
test located in `${BASEPATH}/lmfit2/C/tests`:

.. code-block:: bash
    cd ${BASEPATH}/lmfit2/C/tests
    bash test.bash

This might take a minute, but you shouldn't end up with any errors if
everything is working as it should.


Python Version Installation
---------------------------

The python version depends on the `lmfit` python library. You can install it
using `pip`, like so:

.. code-block:: bash
    pip install lmfit


The python version is not set up to be a python package, instead it is
intened to be used as a script. To "install" it, simply copy both the
`lmfit2.py` and `dmap.py` python files to whereever you want to them to be
to use them.


Test the Python Code
--------------------

Assuming that you have grabbed the source code for lmfit2 using `git`
(see the above section: C Version Installation), navigate to the 
`lmfit2/python` directory and run:

.. code-block:: bash
    python tests/test.py

This will take a minute or two. It shouldn't return any errors if everything
is set up properly.