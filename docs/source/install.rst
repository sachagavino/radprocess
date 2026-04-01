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

Inside radprocess/ run::

    ./setup.sh

Follow the information prompted on the terminal. If everything works fine, RadProcess is now installed and ready to be used by typing anywhere in your terminal::

    radprocess  

Alternatively, you can do::

    source .venv/bin/activate

And then run::

    radprocess

Because radprocess uses specific softwares (RADMC3D, POLARIS, PYMSES), it is necessary that these softwares are preliminary installed in the user's machine 
before using RadProcess.



Running the code
=================

The code is simply run by typing::

    radprocess

RadProcess is embedded in a Gradio interface, which can be launched by clicking on one of the links that appear. If you run radprocess on a local machine, it is recommanded to
clic on the local link. If, on the other hand, you use RadProcess from a server, it is required that you use the public link instead. 
