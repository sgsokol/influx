
.. _install:

.. highlight:: bash


============
Installation
============

To use the software ``influx_si``, you'll need some
dependencies listed bellow. The software was developed on Linux
but can be used both on Linux (or other UNIX, MacOS included) and Windows platforms.
If you are not used to install system wide environments
like R or Python, ask your local computer
support for help. We don't provide support for installation.

.. note:: The code examples here after are given for Unix shell environment.
 On windows, in DOS environment the syntax is often similar and in
 cygwin or Ubuntu environment (Unix tools on Windows) the syntax is identical
 to the Unix's one.


Dependencies
------------

- R-3.4.0 or higher (cf http://www.r-project.org/ or your system packaging solution) + the following packages.
  
  + nnls
  + rmumps (5.2.1-3 or higher)
  + arrApply
  + slam
  + limSolve (optional, needed only for ``--lim`` option)

  To install R modules, as administrator do in R

  .. code-block:: rconsole

     > install.packages(c("nnls", "rmumps", "arrApply", "slam", "limSolve"), dep=TRUE)
 
  (you can adapt the package list according to your needs by removing optional packages)

  If you are not an administrator of your R installation, you can execute the command above in your own session and install necessary packages in your own disk space. Other users will have to do the same in their respective sessions if they want to use ``influx_si``.

- python 2.6 (or higher but not 3.0 or higher) and modules

  + numpy
  + libsbml (optional, needed for ftbl2metxml.py)
- cytoscape is optional (http://www.cytoscape.org).
  It can be used to visualize your networks
  by intermediate of ``ftbl2xgmml.py`` utility.
  You can also map flux values returned by ``influx_si`` on some
  graphical parameter like edge width for visualizing purposes.

Python and R are advised to be in your PATH variable,
in other words, they should be executable from any directory.

.. warning:: As of this writing (September 17, 2014), an R package ``nnls`` distributed in precompiled form on Windows platform, can produce wrong results if a 32 bits version is used on Windows 64 bits. To avoid this, use 64 bit version of R on Windows 64 bits or recompile it by hand. To be sure to use 64 bits version of R, check that the ``Path`` system variable has the R path ending by ``\bin\x64`` and not just by ``\bin``.

.. note:: On some Python distributions (e.g. Anaconda) on Windows platform, association between ``.py`` files and Python interpreter is made in incomplete way: the file is executed but command line arguments are not passed to Python. To correct this, a user with administrator privileges has to edit register base with ``regedit``. The key ``HKEY_CLASSES_ROOT\py_auto_file\shell\open\command`` must be changed from
  
   .. code-block:: text
   
     "<path_to_your_python_dir>\python.exe" "%1"
  
   to
   
   .. code-block:: text
   
     "<path_to_your_python_dir>\python.exe" "%1" %*


   It may happen (depending on your Windows version) that some other keys (related to Python too) have to be modified in similar way.

Compilation dependencies
------------------------

Starting from version 3.0, some critical parts of code are written in C++ which will require a corresponding compiler installed on your system. It is strongly advised to use the same compiler that was used to compile your R software. You can find which one it was by checking the output of the following shell command ::

$ R CMD config CXX

It is likely to be ``g++``. A compilation for a given version of ``influx_si`` will be done automatically only once at the very first execution of ``influx_s.py`` or ``influx_i.py``.

On Linux, all tools necessary for compilation are often available by default. If not, install Linux package (as well as its dependencies) containing ``g++`` compiler (or what ever was used to compile R).

If you are on Windows platform, you have to install RTOOLS software collection available from https://cran.r-project.org/bin/windows/Rtools/
Be sure to pick up a frozen version that corresponds to your R version. This package will contain the necessary C++ compiler.

If you are on MacOS, your have to install Xcode from AppStore. Furthermore, if some of required R packages are not available in binary form for installation, they will be compiled from sources and this can require additional installation of gfortran-4.8 (or higher).

``influx_si`` installation
--------------------------

Unpack the content of ``influx_si-vX.Y.zip`` (where X.Y is the version number)
somewhere on your disk. If you want to make ``influx_si`` available
system wide and install it in a protected directory, you need
administrative privileges. Otherwise, ``influx_si`` will be
available only in your personal session.

Add this new directory to your (or system wide) PATH variable
(if you don't know what does it mean or how to do it,
ask for help from your local computer service).
This step is optional but if you don't do it, you
need to type all the path to ``influx_si`` and their utilities
every time you run it. It can be as cumbersome as ::

$ /home/joe/soft/bio/flux/influx_s-v2.9/influx_s.py mynetwork.ftbl

instead of simple ::

$ influx_s.py mynetwork.ftbl

If you want to make ``influx_si`` available system wide without
modifying the PATH variable, add a symbolic link in a directory
which is already in PATH. For example, as root you can do

:: 

  $ cd /usr/local/bin
  $ ln -s /path/to/dir/of/influx_s/{influx_s.py,influx_i.py,res2ftbl_meas.py,ftbl2cumoAb.py,ftbl2kvh.py,ftbl2netan.py,ftbl2xgmml.py,ff2ftbl.py,ffres2ftbl.py,txt2ftbl.py,ftbl2metxml.py} .

assuming that ``/usr/local/bin`` is already in the PATH.

First compilation
-----------------
To accomplish the installation, you have to run ``influx_s.py`` or ``influx_i.py`` for the first time as a user having write permissions to the installation directory. I.e. if you have installed ``influx_si`` as system administrator you have to make a first run also as a system administrator. This first run will compile a shared library ``mult_bxxc.so`` (a suffix ``.so`` can be different on your platform) needed for further ``influx_si`` executions. An example of a command to run is given in the next session "Test of installation".

If in the future, for any reason (upgrading R version, changing the compiler, ...) you have to recompile the shared library, just remove the file ``mult_bxxc.so`` (or its equivalent if you are not on a Linux platform) and rerun ``influx_si`` on any FTBL file being a user with write permission on installation directory.

********************
Test of installation
********************
Open a shell window and set your current directory to the ``<influx_si_install_dir>/test``.
To run ``influx_s`` you can type ::

 $ influx_s.py e_coli.ftbl

or ::

 $ ../influx_s.py e_coli.ftbl

if it is not in the PATH

or drag-and-drop the icon of ``e_coli.ftbl`` to the icon of ``influx_s.py``.

If everything was correctly installed, you should see in your shell window an
output looking like:

.. code-block:: text

 "../influx_s.py" "e_coli.ftbl"
 code gen: 2016-07-29 12:06:04
 calcul  : 2016-07-29 12:06:04
 end     : 2016-07-29 12:06:08

The meaning of this output is quit simple. First, an R code is generated from FTBL file then it is executed till it ends. Time moments at which these three events occur are reported.

At the very first execution, a compilation of auxiliary file ``mult_bxxc.cpp`` will occur which will modify the output in the following manner

.. code-block:: text

 "../influx_s.py" "e_coli"
 code gen: 2016-04-12 10:45:31
 calcul  : 2016-04-12 10:45:31
 g++ -I/usr/local/src/R-3.2.4/include -DNDEBUG  -I/usr/local/include  -I"/home/local/src/R-3.2.4/library/Rcpp/include" -I"/home/local/src/R-3.2.4/library/RcppArmadillo/include" -I"/home/local/src/R-3.2.4/library/rmumps/include" -I"/home/sokol/insa/sysbio/dev/ftbl2sys"    -fpic  -O2 -mtune=native -ffast-math  -O3 -mtune=native -std=c++11 -c mult_bxxc.cpp -o mult_bxxc.o
 g++ -shared -L/usr/local/src/R-3.2.4/lib -L/usr/local/lib64 -o sourceCpp_1.so mult_bxxc.o -L/usr/local/src/R-3.2.4/lib -lRlapack -L/usr/local/src/R-3.2.4/lib -lRblas -lgfortran -lm -lquadmath /home/local/src/R-3.2.4/library/rmumps/libs/rmumps.so -L/usr/local/src/R-3.2.4/lib -lRlapack -L/usr/local/src/R-3.2.4/lib -lRblas -lgfortran -lm -lquadmath -L/usr/local/src/R-3.2.4/lib -lR
 end     : 2016-04-12 10:45:44

On your system, the compilation commands and paths can differ from this example. That's normal.

The calculation result will be written in ``e_coli_res.kvh``.
It should be almost identical to the same file in ``ok/`` subdirectory.
On Unix you can do ::

$ diff e_coli_res.kvh ok/e_coli_res.kvh

to see if there is any difference. Some small differences in numerical
values can be ok. They might come from variations in versions of R and
underlying numerical libraries (BLAS, LAPACK and so on).

If something went wrong, check the error messages in ``e_coli.err``,
interpret them, try to figure out why the errors occurred and correct them.

In high throughput context, you can find useful to run ``influx_si`` in parallel on many FTBL files. It can be done just by providing more than one FTBL file in argument. For example, with two of FTBLs provided with the package you can run: ::
 
 $ ../influx_s.py e_coli.ftbl e_coli_growth.ftbl
 

In this case, the output looks sightly different than in one by one run:

.. code-block:: text

 "../influx_s.py" "e_coli.ftbl" "e_coli_growth.ftbl"
 e_coli: code gen: 2016-07-29 12:13:32
 e_coli_growth: code gen: 2016-07-29 12:13:32
 //calcul: 2016-07-29 12:13:32
 //end   : 2016-07-29 12:13:36
 
The time moments for code generation is preceded by a short version of FTBL file names. The symbol ``//`` means parallel proceeding. Parallel calculations are launched after all files are proceeded for the code generation.

It is the operating system that dispatches and equilibrates the charge
among available CPUs and cores, not ``influx_si`` who simply launches these processes.

For a quick test of ``influx_i``, you can run in the same directory ::

$ ../influx_i.py e_coli_i

Normal output looks like

.. code-block:: text

 "../influx_i.py" "e_coli_i"
 code gen: 2016-04-12 10:43:10
 calcul  : 2016-04-12 10:43:10
 end     : 2016-04-12 10:43:35

Calculation results are written in ``e_coli_i_res.kvh`` and they can be compared with the same file in the ``ok/`` sub-directory. You can also visually check a generated graphic file ``e_coli_i.pdf`` to see if all simulated label kinetics based on estimated fluxes and metabolite concentrations are close to experimental data.

For a quick start guide, launch ::

$ influx_s.py --help

or ::

$ influx_i.py --help

depending on what context you want to treat: stationary or instationary labeling.

These commands show all available options with a brief description.
For more detailed documentation read :doc:`User's manual <manual>`.
