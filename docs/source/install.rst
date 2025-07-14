Installation of radprocess
************

How to obtain radprocess
=================

radprocess can be obtained by cloning its repository (soon using pip). This is also the best way to keep up-to-date.
From a terminal, go to a directory where you want to install the code, and type:: 


    git clone https://github.com/sachagavino/radprocess.git


This will create a folder ``radprocess/``, which contains the full git repository. You can now access the package::


    cd radprocess/


To make sure you always use the latest version, you can type:: 


    git pull



Requirements and environment
=================

Because radprocess (for now) uses specific library versions, it is strongly recommended to use a dedicated virtual environment to avoid conflicts with other packages. You can use ``conda`` to create a virtual environment. 
The easiest way to do this is to use the provided ``environment.yml`` file. The name of the environment is ``radprocess`` by default, but you can change it in the ``environment.yml`` file before creating the environment.
From the terminal, type::

    conda env create -f environment.yml

Verify that the environment is created properly::

    conda env list

If the name of the environment appears in the list, it means the environment has been created successfully. 

You can now activate the new environment with::

    conda activate radprocess





Running the code
=================


