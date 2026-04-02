Installation of radprocess
**************************

How to obtain radprocess
========================

The radprocess package can be obtained by cloning its GitHub repository
(pip installation will be available in the future). Cloning the repository
is currently the recommended method, as it allows you to easily keep your
installation up to date.

From a terminal, navigate to the directory where you wish to install the code
and run::

    git clone https://github.com/sachagavino/radprocess.git

This will create a folder named ``radprocess/`` containing the full repository.
You can then move into the directory::

    cd radprocess/

To update your local copy to the latest version, simply run::

    git pull


Requirements and environment
============================

From within the ``radprocess/`` directory, run::

    ./setup.sh

Follow the instructions displayed in the terminal. If the setup completes
successfully, RadProcess will be installed and ready to use.

You can then launch the code from anywhere in your terminal with::

    radprocess

Alternatively, you may activate the virtual environment manually::

    source .venv/bin/activate

and then run::

    radprocess

RadProcess relies on several external software packages, including
RADMC-3D, POLARIS, and PyMSES (Python>=3). These must be installed and properly
configured on your system prior to using RadProcess. 


Running the code
================

RadProcess is launched simply by typing::

    radprocess

This will start a Gradio-based graphical interface. One or more URLs will be
displayed in the terminal.

- On a local machine, it is recommended to open the **local URL**.
- When running on a remote server, you should use the **public (share) URL**
  if available. If the share link is not available (e.g., due to network
  restrictions), you can access the interface via SSH tunneling.

In that case, run the following command from your local machine::

    ssh -L 7860:localhost:7860 user@server

and then open the following address in your browser::

    http://localhost:7860
