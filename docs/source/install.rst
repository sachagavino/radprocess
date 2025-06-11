Installation of radprocess
************

How to obtain radprocess
=================

radprocess can be obtained by cloning its repository. This is also the best way to keep up-to-date.
From a terminal, go to a directory where you want to install the code, and type:: 


    git clone https://github.com/sachagavino/radprocess.git


This will create a folder ``radprocess/``, which contains the full git repository. You can now access the package::


    cd radprocess/


To make sure you always use the latest version, you can type:: 


    git pull



Requirements
=================

The following softwares are required:

#. ``make``

    The GNU Make tool is required to compile the software. In principle, it should already be pre-installed on your machine.
    In case it is not, you can type ``sudo apt-get install make`` if you are working on Linux, or ``brew install make`` regardless of your OS.


#. ``Fortran-90 compiler``

    You need a Fortran-90 compiler. The software has been tested with the ``gfortran`` compiler only, but there is no reason it should not be working with the others. Please, make sure you have the latest version
    of your compiler (gfortran version > 10.3.0). 


Running the code
=================


