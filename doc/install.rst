
.. _install:


============
Installation
============

To install the software ``influx_s``, you'll need some
dependencies listed bellow. The software was developed on Linux
but can be used both on Linux (or other UNIX, MacOS included but
not yet tested) and Windows platforms.
If you are not used to install system wide environments
like R or Python, ask some help from your local computer
service. We don't provide support for installation.

.. note:: The code examples here after are given for Unix shell environment.
 On windows, in DOS environment the syntax is often similar and in
 cygwin environment (Unix tools on Windows) the syntax is identical
 to the Unix's one.


Dependencies
------------

- R-2.12 (or higher, cf http://www.r-project.org/ or your system packaging solution) + the following packages
  
  + bitops
  + fUtilities
  + Matrix
  + nnls
  + optparse
  + multicore (if you plan to use Monte-Carlo on Unix platform. On Windows this option is not working well yet)
  + Rtools (on Windows, you can get it from http://www.murdoch-sutherland.com/Rtools/)
    or any other mean to make the command ``R CMD SHLIB --clean`` working on fortran files.
To install R modules, as administrator do in R::

 > install.packages(c("fUtilities", "bitops", "Matrix", "nnls", "optparse"), dep=T)
 > install.packages("multicore",,"http://www.rforge.net/")

- python 2.6+ (not 3.0 or higher) and module
  
  + numpy
- cytoscape is optional (http://www.cytoscape.org).
  It can be used to visualize your networks
  by intermediate of ftbl2cytoscape utility.
All these components are advised to be in your PATH variable,
in other words, they should be executable from any directory.

influx_s installation
---------------------
Unpack the content of ``influx_s-vX.Y.tgz`` or .zip (where X.Y is the version number)
somewhere on your disk. If you want to make ``influx_s`` be available
system wide and install it in a protected directory, you need
administrative privileges for this. Otherwise, ``influx_s`` will be
available only in your personal session.

Add this new directory to you (or system wide) PATH variable
(if you don't know what does it mean or how to do it,
ask help from your computer service).
This step is optional but if you don't do it, you
need to type all the path to ``influx_s`` every time that you run
it. It can be as cumbersome as ::

$ /home/joe/soft/bio/flux/influx_s/influx_s.py mynetwork.ftbl

instead of simple ::

$ influx_s.py mynetwork.ftbl

If you want to make influx_s available system wide without
modifying the PATH variable, add a symbolic link in a directory
which is already in PATH. For example, as root you can do ::

$ cd /usr/local/bin
$ ln -s /path/to/dir/of/influx_s/influx_s.py .

assuming that /usr/local/bin is in the PATH.

**********
Test of installation
**********
Open a shell window and set your current directory
to the ``<influx_s_install_dir>/test``.
To run influx_s you can type ::

$ influx_s.py e_coli_lcms.ftbl

or drag-and-drop the icon of ``e_coli_lcms.ftbl`` to the icon of ``influx_s.py``

If everything was correctly installed, you should see in your shell window an
output looking like: ::

 code gen: 2011-06-28 14:54:22
 calcul  : 2011-06-28 14:54:27
 end     : 2011-06-28 14:54:44

The result file will be in ``e_coli_lcms_res.kvh``.
It should be identical to the same file in ``ok/`` subdirectory.
On Unix you can do ::

$ diff e_coli_lcms_res.kvh ok/e_coli_lcms_res.kvh

to see if there is any difference. Some small differences in numerical values can be ok. They might come from variations in versions of R and depending numerical libraries (BLAS, LAPACK and so on).

If something get wrong, check the error messages in ``e_coli_lcms.err``,
interpret them and try to correct the underlying errors.

For quick start guide, launch ::

$ influx_s.py --help

it shows all available option with brief description.
Or read :doc:`User's manual <manual>`  for more details.
