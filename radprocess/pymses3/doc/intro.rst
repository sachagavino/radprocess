####################################
PyMSES : Python modules for RAMSES_
####################################

************
Introduction
************

PyMSES is a set of Python modules originally written for the RAMSES_ astrophysical fluid dynamics AMR code.

Its purpose :

  * provide a clean and easy way of **getting the data** out of RAMSES_ simulation outputs.
  * help you analyse/manipulate very large simulations transparently, without worrying more than needed about domain
  decompositions, memory management, etc.,
  * interface with a lot of powerful Python libraries (Numpy_/Scipy_, Matplotlib_, PIL_, HDF5_/PyTables_) and existing
  code (like your own).
  * be a post-processing toolbox for your own data analysis.

.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://www.scipy.org/
.. _PIL: http://www.pythonware.com/products/pil/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PyTables: http://www.pytables.org/

What PyMSES is NOT
==================

It is **not an interactive environment** by itself, but :

  * it provides modules which can be used interactively, for example with IPython_.
  * it also provides an :doc:`AMRViewer <ug_visu_AMRViewer>` graphical user intergace (GUI) module that allows you to
  get a quick and easy look at your AMR data.


*************
Downloads
*************

* :doc:`downloads`


*************
Documentation
*************

* :doc:`index`


********
Contacts
********

    :Authors:       Damien CHAPON, Marc LABADENS, Thomas GUILLET, Olivier IFFRIG
    :Contact:       damien.chapon@cea.fr
    :Organization:  `IRFU, CEA/Saclay <http://irfu.cea.fr/>`_
    :Address:       Gif-Sur-Yvette, F91190, France
    :Date:          |today|

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _RAMSES: http://www.itp.uzh.ch/~teyssier/Site/RAMSES.html
.. _IPython: http://ipython.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _mayavi: http://code.enthought.com/projects/mayavi/

