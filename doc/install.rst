
.. _install:


============
Installation
============

To install the software ``influx_s``, you'll need some
dependencies listed bellow. The software was developed on Linux
but can be used both on Linux (or other UNIX, MacOS included but
not yet tested) and Windows platforms.
If you are not used to install system wide environments
like R or Python, ask your local computer
support for help. We don't provide support for installation.

.. note:: The code examples here after are given for Unix shell environment.
 On windows, in DOS environment the syntax is often similar and in
 cygwin environment (Unix tools on Windows) the syntax is identical
 to the Unix's one.


Dependencies
------------

- R-2.15 (or higher, cf http://www.r-project.org/ or your system packaging solution) + the following packages
  
  + bitops
  + nnls

To install R modules, as administrator do in R::

 > install.packages(c("bitops", "nnls"), dep=T)

If you are not an administrator of your R installation, you can execute the command above in your own session and install necessary packages in your own disk space. Other users will have to do the same install in their respective sessions if they want to use ``influx_s``.

- python 2.6+ (not 3.0 or higher) and module
  
  + numpy
- cytoscape is optional (http://www.cytoscape.org).
  It can be used to visualize your networks
  by intermediate of ftbl2cytoscape utility.
  You can also map flux values returned by ``influx_s`` on some
  graphical parameter like edge width for visualizing purposes.

All these components are advised to be in your PATH variable,
in other words, they should be executable from any directory.

influx_s installation
---------------------
Unpack the content of ``influx_s-vX.Y.zip`` (where X.Y is the version number)
somewhere on your disk. If you want to make ``influx_s`` available
system wide and install it in a protected directory, you need
administrative privileges. Otherwise, ``influx_s`` will be
available only in your personal session.

Add this new directory to your (or system wide) PATH variable
(if you don't know what does it mean or how to do it,
ask for help from your computer service).
This step is optional but if you don't do it, you
need to type all the path to ``influx_s`` every time you run
it. It can be as cumbersome as ::

$ /home/joe/soft/bio/flux/influx_s-v2.2/influx_s.py mynetwork.ftbl

instead of simple ::

$ influx_s.py mynetwork.ftbl

If you want to make ``influx_s`` available system wide without
modifying the PATH variable, add a symbolic link in a directory
which is already in PATH. For example, as root you can do ::

$ cd /usr/local/bin
$ ln -s /path/to/dir/of/influx_s/influx_s.py .

assuming that ``/usr/local/bin`` is already in the PATH.

********************
Test of installation
********************
Open a shell window and set your current directory
to the ``<influx_s_install_dir>/test``.
To run ``influx_s`` you can type ::

 $ influx_s.py e_coli.ftbl

or ::

 $ ../influx_s.py e_coli.ftbl

if it is not in the PATH

or drag-and-drop the icon of ``e_coli.ftbl`` to the icon of ``influx_s.py``

If everything was correctly installed, you should see in your shell window an
output looking like: ::

 "../influx_s.py" "e_coli.ftbl"
 code gen: 2013-02-15 16:42:37
 calcul  : 2013-02-15 16:42:44
 end     : 2013-02-15 16:43:06

The result file will be in ``e_coli_res.kvh``.
It should be almost identical to the same file in ``ok/`` subdirectory.
On Unix you can do ::

$ diff e_coli_res.kvh ok/e_coli_res.kvh

to see if there is any difference. Some small differences in numerical
values can be ok. They might come from variations in versions of R and
underlying numerical libraries (BLAS, LAPACK and so on).

If something get wrong, check the error messages in ``e_coli.err``,
interpret them and try figure out and correct corresponding errors.

In high throughput context, you can find usefull to run ``influx_s`` in parallel
on many FTBL files. It can be done just by providing more than one FTBL file in argument. For example, with two FTBL provided with the package you can run: ::

 $ ../influx_s.py e_coli.ftbl e_coli_growth.ftbl

In this case, the output looks sightly different than in one by one run: ::

 "../influx_s.py" "e_coli" "e_coli_growth"
 e_coli: code gen: 2013-02-22 18:23:10
 e_coli_growth: code gen: 2013-02-22 18:23:10
 e_coli_growth: calcul  : 2013-02-22 18:23:13
 e_coli: calcul  : 2013-02-22 18:23:18
 e_coli: end     : 2013-02-22 18:23:42
 e_coli_growth: end     : 2013-02-22 18:23:48

The time moments for code generation, calculation and end of program
are preceded by a short version of FTBL file names. As all files are
proceeded in parallel, there is no a predefined order in which all these
records can appear in the log console.

It is the operating system that dispatches and equilibrates the charge
among available CPUs and cores, not ``influx_s`` who simply launches these processes.

For quick start guide, launch ::

$ influx_s.py --help

it shows all available option with a brief description.
For more detailed documentation read :doc:`User's manual <manual>`.
