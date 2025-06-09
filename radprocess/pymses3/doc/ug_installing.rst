Installing PyMSES
=================

Requirements
------------

PyMSES has some :ref:`core-dep` plus some :ref:`recommended-dep` you might need to install to enable all PyMSES features.

**The developpment team** strongly recommends the user to install the EPD_ (Enthought Python Distribution) which wraps all these dependencies into a single, highly-portable package.

.. _EPD: http://www.enthought.com/products/epd.php

.. _core-dep:

Core dependencies
.................

These packages are mandatory to use the basic functionality of PyMSES:

  * a gcc-compatible C compiler,
  * GNU make and friends, 
  * Python_, version 2.x (*not* 3.x), *including developement
    headers* (Python.h and such), python 2.6.x or later is recommended to use some multiprocessing speed up.
  * Python packages:

    - numpy_, version 1.2 or later
    - scipy_

  * iPython_ is not strictly required, but it makes the interactive experience so much better you will certainly want to install it.

.. _Python: http://www.python.org/
.. _numpy: http://numpy.scipy.org/
.. _scipy: http://www.scipy.org/
.. _iPython: http://ipython.scipy.org/

.. _recommended-dep:

Recommended dependencies
........................

Those packages are recommended for general use (plotting, easy input and output, image processing, GUI, ...). Some PyMSES features may be unavailable without them:

  * matplotlib_ for plotting
  * the Python Imaging Library (PIL_) for Image processing
  * HDF5_ and PyTables_ for Python HDF5 support
  * wxPython_ for the AMRViewer GUI
  * mpi4py_ if you want to use the MPI library on a large parallel system.

.. _wxPython: http://www.wxpython.org/
.. _PIL: http://www.pythonware.com/products/pil/
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PyTables: http://www.pytables.org/
.. _mpi4py: http://code.google.com/p/mpi4py/

Delevoper dependencies
......................

You will need this if you intend to work on the source code, or if you obtained PyMSES for an unpackaged version (i.e. a tarball of the mercurial repository, or ``hg clone``)

  * mercurial_ for source code management
  * Cython_
  * sphinx_ for generating the documentation

.. _mercurial: http://mercurial.selenic.com
.. _Cython: http://www.cython.org/
.. _sphinx: http://sphinx.pocoo.org/

Installation instructions
-------------------------

For now, the easiest way to install PyMSES from a tarball is:

  #. Extract the tarball into a directory, say ~/codes/pymses
  #. Run make in the directory
  #. Add the make directory to your PYTHONPATH
  #. Optional : Add the pymses/bin to your PATH, to quickly start the GUI with the amrviewer command or to launch basic scripts.

For example, using the bash shell:

.. code-block:: bash

    $ cd  ~/codes
    $ tar -xvfz pymses-3.0.0.tar.gz
    $ cd pymses_3.0.0
    $ make
    $ export PYTHONPATH=~/codes/pymses_3.0.0:$PYTHONPATH
    $ export PATH=$PATH:~/codes/pymses_3.0.0/bin

Note that you will need to place the export statements in your :file:`~/.bashrc` or equivalent to set your ``PYTHONPATH`` and ``PATH`` for all future shell sessions.


