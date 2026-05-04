.. _chap-install:

Installation of radprocess
**************************

Requirements
============

- **Python 3.12** (3.11 may work, but 3.13+ is not yet supported due to pymses)
- A working installation of `pymses <https://irfu.cea.fr/Projets/PYMSES/>`_
  (provided internally, see below)
- Optional: `POLARIS <https://portia.astrophysik.uni-kiel.de/polaris/>`_ and
  `RADMC-3D <https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_
  executables, needed only for running Steps 4, 6, and 8 of the pipeline


How to obtain radprocess
========================

Clone the GitHub repository, create a virtual environment, and install in
editable mode::

    git clone https://github.com/sachagavino/radprocess.git
    cd radprocess
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e .

This creates a virtual environment (``.venv/``) inside the ``radprocess/``
directory, activates it, and installs radprocess with all its Python
dependencies (numpy, scipy, zarr, matplotlib, astropy, etc.).

.. note::

   Always use ``python -m pip`` instead of bare ``pip``. On some systems,
   ``pip`` may point to a different Python environment than the one activated
   by the venv, which leads to packages being installed in the wrong place.

If you also want the Gradio web interface::

    python -m pip install -e ".[gui]"

Every time you open a new terminal and want to use radprocess, you need to
activate the environment first::

    cd /path/to/radprocess
    source .venv/bin/activate

To update your local copy later::

    git pull
    python -m pip install -e .


Alternative: using uv
---------------------

If you have `uv <https://docs.astral.sh/uv/>`_ installed::

    git clone https://github.com/sachagavino/radprocess.git
    cd radprocess
    uv venv --python 3.12
    source .venv/bin/activate
    uv pip install -e .


Installing pymses
=================

pymses is not available on PyPI and must be installed from the local source
directory. The up-to-date version lives on the server (typically in the shared
software folder). To install it, activate the radprocess environment first,
then install pymses into it::

    # Activate the radprocess environment
    source /path/to/radprocess/.venv/bin/activate

    # Go to the pymses source directory on the server
    cd /path/to/pymses

    # Install pymses (--no-build-isolation ensures it uses the numpy
    # and Cython already installed in the environment)
    python -m pip install -e . --no-build-isolation

This compiles the Cython extensions against the correct numpy version and
registers pymses so that Python can find it.

.. warning::

   Do **not** use ``$PYTHONPATH`` or ``$PATH`` as an alternative.

   - ``$PATH`` only affects executables, not Python imports.
   - ``$PYTHONPATH`` bypasses pip, can cause version conflicts across
     projects, and may lead to crashes if pymses was compiled against a
     different numpy.


Using radprocess from a notebook
=================================

Make sure the notebook kernel uses the same Python environment where radprocess
is installed. You can verify from a notebook cell:

.. code-block:: python

    import sys
    print(sys.executable)

and compare with the output of ``which python`` in the terminal where you ran
``python -m pip install -e .``.


Registering the environment as a Jupyter kernel
-----------------------------------------------

If your notebook uses a different default kernel, register the radprocess
environment explicitly::

    # From the terminal, with the radprocess venv activated:
    python -m pip install ipykernel
    python -m ipykernel install --user --name radprocess --display-name "Python (radprocess)"

Then select **Python (radprocess)** as the kernel in your notebook.


Running the Gradio interface
============================

If installed with the ``[gui]`` extra, radprocess can be launched as a
Gradio-based graphical interface::

    radprocess

One or more URLs will be displayed in the terminal.

- On a local machine, open the **local URL**.
- On a remote server, use the **public (share) URL** if available, or set up
  SSH tunneling::

      ssh -L 7860:localhost:7860 user@server

  then open ``http://localhost:7860`` in your browser.


Verifying the installation
==========================

.. code-block:: python

    from radprocess.pipeline.Pipeline import Pipeline

    pipe = Pipeline()
    pipe.configparams.ramsesoutput.ramses_output_dir = "/path/to/output_00940/"
    pipe.set_working_dir("/path/to/working_directory")

    # Check that RAMSES output is readable:
    print(pipe.read_hydro_descriptor())
    print(pipe.read_sink_info())


Troubleshooting
===============

**"No module named pymses"**
    pymses must be installed separately into the radprocess environment.
    See the section above.

**zarr import errors or "BloscCodec not found"**
    radprocess requires zarr v3. Check your installed version::

        python -c "import zarr; print(zarr.__version__)"

    If it shows 2.x, upgrade with ``python -m pip install "zarr>=3.0"``.

**"No module named gradio"**
    Gradio is optional. Install it with ``python -m pip install -e ".[gui]"``
    or ``python -m pip install gradio``. It is only needed for the web
    interface, not for notebook or script usage.

**numpy/Cython build errors when installing pymses**
    pymses requires Cython at build time. Make sure Cython is installed
    (``python -m pip install cython``) before installing pymses. If you still
    get errors about deprecated numpy C API, the internal pymses fork should
    handle this.

**pip points to the wrong environment**
    If ``which pip`` shows a path outside ``.venv/``, always use
    ``python -m pip`` instead. This guarantees you are using the pip that
    matches the active Python interpreter.