
.. _quick:

===========
Quick Start
===========

A basic work-flow with ``influx_s`` is composed of the following steps:

1. Create a FTBL file describing your metabolic reactions, carbon transitions, experimental data and some options. Let call an example file ``mynetwork.ftbl``. The FTBL file must follow syntax rules elaborated for `13CFlux <https://www.13cflux.net/>`_ software. The FTBL file is a plain text file. The syntax rules will be more or less obvious for someone working on metabolism biochemistry. So, to go quickly, you can inspire from an example file ``test/e_coli.ftbl`` distributed with the ``influx_s`` software.
 .. note:: Staring from the version 2.5, ``NA`` values (as "Non Available") are admitted as measurements values where appropriate. The difference with FTBL where they are simplly omitted is that NA measurments are simulated and are present in the vectors ``simulated unscaled labeling measurements`` and ``simulated scaled labeling measurements`` in the result kvh file.

2. Set your current directory to the directory of ``mynetwork.ftbl`` and run::

 ``$ influx_s.py mynetwork``

 or::

  $ /path/to/install/dir/of/influx_s/influx_s.py mynetwork

 Note that the suffix ``.ftbl`` is optional.

 The ``influx_s`` run will produce the following files in the same directory that ``mynetwok.ftbl``

 ``mynetwork.log``
   containing the run-time output from various scripts, in particular,
   it contains a report on convergence history during the fitting process.
   It can be helpfull for identifying potential problems but if everything
   is going well, the user does not have to examine the content of this file;
 ``mynetwork.err``
  containing the warning and error messages.
  Normally, this file should be empty (0 byte size);
 ``mynetwork_res.kvh``
  containing all of the results. `KVH format <http://serguei.sokol.free.fr/kvh-format/>`_ is a
  lightweight plain text format for hierarchically structured data. It can be seen in a text editor
  or in a spreadsheet software as its fields are tab separated. It can also be processed by user's
  custom software for post-processing, graphics output and alike. If ``influx_s``
  is run on a series of starting points there will be generated a result
  file ``mynetwork_res.kvh`` containing common information to all starting points
  but also one kvh file by starting point, e.g. ``mynetwork_res.V1.kvh``,
  ``mynetwork_res.V2.kvh`` and so on;
 ``mynetwork.pres.txt``
  containing a matrix of fitted parameters and final cost values. Each column
  corresponds to a particular starting point if run with ``--fseries`` and /or
  ``--iseries`` options. If ``influx_s`` was run without these options, the file
  will contain only one column corresponding to the starting point defined
  in the ``mynetwork.ftbl`` file.
  
 ``edge.netflux.mynetwok``, ``edge.xchflux.mynetwok``, ``node.log2pool.mynetwork``
  as the middle name of this files suggest, they can be used to map the corresponding
  values on the network graph in the `cytoscape <http://www.cytoscape.org>`_ software.

  .. note:: All these files are silently overwritten if already exist.
   So take care to copy your results elsewhere if you want to protect them
   from overwriting.

 .. note:: It can be helpful to do some "dry runs" by executing ::

   $ influx_s.py --noopt mynetwork
   
   before collecting actual data measurement to see if intended measurements will be sufficient to well define all fluxes or at least the fluxes of interest. It is possible to do because the measurement values in the FTBL file does not matter for flux SD calculation when ``--noopt`` option is used. So it can be used any values even NA at this moment. In the contrary, ``dev`` values set in the FTBL file, must be realistic. It is generally not a problem as they express measurements errors and are more or less known for a given measurement chain.
   
   It is worthwhile to stress that a "dry run" is done for some presumed free fluxe values and if they reveal to be very different from actual flux values, it can happen that a network considered as well defined at moment of "dry run" turned into a badly defined network with actual measurement data and corresponding estimated fluxes. So it is important to do his best to guess the most realistic free fluxes for "dry runs".

3. See warning and error messages in ``mynetwork.err`` if any. Correct what has to be corrected and retry p. 2

4. Extract and use the numerical results from the ``mynetwork_res.kvh`` file.

5. Optionally, visualize net fluxes (or exchange fluxes or :math:`\log_2(M)`) in cytoscape using ``edge.netflux.mynetwok``, ``edge.xchflux.mynetwok`` or ``node.log2pool.mynetwork``.
