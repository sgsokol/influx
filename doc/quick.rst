
.. _quick:

===========
Quick Start
===========

A basic work-flow with ``influx_s`` is composed of the following steps:

1. Create a FTBL file describing your metabolic reactions, carbon transitions, experimental data and some options. Let call an example file ``mynetwork.ftbl``. The FTBL file must follow syntax rules elaborated for `13CFlux <https://www.13cflux.net/>`_ software. The FTBL file is a plain text file. The syntax rules will be more or less obvious for someone working on metabolism biochemistry. So, to go quickly, you can inspire from an example file ``test/e_coli_lcms.ftbl`` distributed with the ``influx_s`` software.

2. Set your current directory to the directory of ``mynetwork.ftbl`` and run:
 ::

  $ influx_s.py mynetwork

 or::

  $ /path/to/install/dir/of/influx_s/influx_s.py mynetwork

 Note that the suffix ``.ftbl`` is optional.

 The ``influx_s`` run will produce the following files in the same directory that ``mynetwok.ftbl``

 ``mynetwork.log``
   containing the run-time output from various scripts, in particular, it contains a report on convergence history during the fitting process. It can be helpfull for identifying potential problems but if everything is going well, the user does not have to examine the content of this file.
 ``mynetwork.err``
  containing the warning and error messages. Normally, this file should be empty (0 byte size).
 ``mynetwork_res.kvh``
  containing all of the results. `KVH format <http://serguei.sokol.free.fr/kvh-format/>`_ is a
  lightweight plain text format for hierarchically structured data. It can be seen in a text editor
  or in a spreadsheet software as its fields are tab separated. It can also be processed by user's
  custom software for post-processing, graphics output and alike.

  .. note:: All these files are silently overwritten if already exist. So take care to copy your results  elsewhere if you want to protect them from overwriting.

3. See warning and error messages in ``mynetwork.err`` if any. Correct them and retry p. 2

4. Extract and use the numerical results from the ``mynetwork_res.kvh`` file

