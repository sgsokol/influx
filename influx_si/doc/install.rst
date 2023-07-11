
.. _install:

.. highlight:: bash

.. _Python: https://www.python.org/

.. _R: https://www.r-project.org/

.. _Anaconda: https://www.anaconda.com/

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html

============
Installation
============

The software was developed on Linux
but can be used both on Linux (or other UNIX, MacOS included) and Windows platforms.

.. note:: The code examples here after are given for Unix shell environment.
 On Windows, in PowerShell or DOS environment the syntax is often similar and in cygwin or Ubuntu environment (Unix tools on Windows) the syntax is identical to the Unix's one.

.. note:: In command examples to run, we use script names with extension `.py`. However, starting from version 5.0.3, this extension can be omitted as all Python_ scripts are doubled with executable files without '.py'. For example, commands: ::

 $ influx_s.py --prefix e_coli
 
 and ::

 $ influx_s --prefix e_coli
 
 are now equivalent. Even if it works on all platforms, it can be particularly useful for Windows, where supplementary effort can be required to associate ``.py`` file with Python interpreter. Using executable programs (i.e. without ``.py`` extension) makes this extra configuration step no more mandatory.

Installation with ``conda``
---------------------------
If you have Anaconda or Miniconda installed on your system, installation of ``influx_si`` resumes to: ::

  $ conda install influx_si -c conda-forge -c bioconda
  
It installs ``influx_si`` itself as well as all needed dependencies both in Python and in R.
  
Installation with ``pip``
-------------------------
If you don't have any version of ``conda`` (neither ``miniconda`` nor ``Anaconda``) but do have a Python and R installed on your system, you can install ``influx_si`` with the following procedure.

You need a python tool called ``pip`` which manages pure python packages. If it is not present on your system, you'll have to install it `first <https://pip.pypa.io/en/stable/installing/>`_ to continue with this method. If you have multiple Python versions installed on your system (e.g. Python2 and Python3) you'll have to use ``pip3`` to install the software in the Python3 universe.

The first step will install only Python part of ``influx_si``: ::

  $ pip3 install influx_si
  
or ::

  $ pip3 install --user influx_si
  
if you wish to install ``influx_si`` not system-wide but only in your own userspace.

To use the software ``influx_si``, you'll need some R dependencies listed below. You can try to install them by: ::

  $ influx_s.py --install_rdep

If this procedure fails, you'll have to solve the underlying problem identified from its error messages and rerun the command again.

R_ dependencies
---------------

As of influx_si version 5.0, user has not to install R dependencies manually from an R session. So they are listed here just for information.

- R-3.4.0 or higher (cf http://www.r-project.org/ or your system packaging solution) + the following packages.
  
  + nnls
  + rmumps (5.2.1-6 or higher)
  + arrApply
  + slam
  + limSolve (optional, needed only for ``--lim`` option)
  + multbxxc
  
.. warning:: As of this writing (September 17, 2014), an R package ``nnls`` distributed in precompiled form on Windows platform, can produce wrong results if a 32 bits version is used on Windows 64 bits. To avoid this, use 64 bit version of R on Windows 64 bits or recompile it by hand. To be sure to use 64 bits version of R, check that the ``Path`` system variable has the R path ending by ``\bin\x64`` and not just by ``\bin``.


Python_ dependencies
--------------------

As of influx_si version 5.0, user has not to install Python dependencies manually. So they are listed here just for information.

- python 3.6 (or higher) and modules

  + scipy
  + libsbml (optional, needed for ``ftbl2metxml.py``)

********************
Test of installation
********************

Open a shell window and get to your working directory.
Copy the distributed test directory to the current directory by running ::

 $ influx_s.py --copy_test
 
then you can get in the newly created directory ``test`` and run some tests:

   .. code-block:: shell
 
     $ cd test/mtf
     $ influx_s.py --prefix e_coli

If everything was correctly installed, you should see in your shell window an output looking like:

.. code-block:: text

 "/home-local/sokol/.local/bin/influx_s" "--prefix" "e_coli"
 code gen: 2022-05-25 12:10:53
 calcul  : 2022-05-25 12:10:53
 end     : 2022-05-25 12:10:55

The meaning of this output is quit simple. First, an R code is generated from input MTF files (cf. :ref:`MTF format <mtf>` for more details) then it is executed till it ends. Time moments at which these three events occur are reported.

The calculation result will be written in directory ``e_coli_res``.
It should be almost identical to the same directory in ``ok/mtf`` subdirectory.
On Unix you can do for example ::

$ diff e_coli_res/e_coli.tvar.sim ../ok/mtf/e_coli_res/e_coli.tvar.sim

to see if there is any difference in estimated fluxes. Some small differences in numerical values can be ok. They might come from variations in versions of R and underlying numerical libraries (BLAS, LAPACK and so on).

If something went wrong, check the error messages in ``e_coli.err``,
interpret them, try to figure out why the errors occurred and correct them.

In high throughput context, you can find it useful to run ``influx_si`` in parallel on many independent MTF sets. It can be done by providing more than one ``--prefix`` options. For example, with two of cases provided with the package you can run: ::
 
 $ influx_s.py --prefix e_coli --prefix e_coli_growth
 

In this case, the output looks sightly different than in one by one run:

.. code-block:: text

  e_coli_growth: code gen: 2022-05-25 14:44:56
  e_coli: code gen: 2022-05-25 14:44:56
  //calcul: 2022-05-25 14:44:57
  //end   : 2022-05-25 14:44:58
 
The time moments for code generation is preceded by a short version of file names. The symbol ``//`` means parallel proceeding. Parallel calculations are launched after all files are proceeded for the code generation.

It is the operating system that dispatches and equilibrates the charge
among available CPUs and cores, not ``influx_si`` who simply launches these processes.

One of the main interest of MTF format introduced in v6.0 is an ability to multiplex constant and variable parts of information describing a set of experiments. In this case, many calculations can run in parallel on inter-dependent input files, cf. ``.vmtf`` description in :ref:`MTF format <mtf>`.

For a quick test of ``influx_i``, you can run in the same directory: ::

  $ influx_i.py --prefix e_coli_i

Normal output looks like

.. code-block:: text

  code gen: 2022-05-25 14:50:51
  calcul  : 2022-05-25 14:50:52
  end     : 2022-05-25 14:51:02

Calculation results are written in ``e_coli_i_res`` directory and they can be compared with the sames files in the ``../ok/mtf/e_coli_i_re`` sub-directory. You can also visually check a generated graphic file ``e_coli_i_res/e_coli_i.pdf`` to see if all simulated label kinetics based on estimated fluxes and metabolite concentrations are close to experimental data.

*****************************
Installation of documentation
*****************************

``influx_si`` is distributed with its documentation. To get it locally accessible from your personal disk space, you can run: ::

 $ influx_s.py --copy_doc

It will create a subdirectory ``doc`` in the current directory. This subdirectory contains ``influx_si.pdf``, all-in-one documentation file but also an ``html`` subdirectory with the documentation viewable in your preferred browser.

The documentation is also available `on-line <https://influx-si.readthedocs.io/>`_.

For a quick reminder of available options, launch ::

$ influx_s.py --help

or ::

$ influx_i.py --help

depending on what context you want to treat: stationary or instationary labeling.

For more detailed documentation, read :doc:`User's manual <manual>`.
