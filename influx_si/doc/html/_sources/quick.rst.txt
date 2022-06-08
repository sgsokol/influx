
.. _quick:

===========
Quick Start
===========

A basic work-flow with ``influx_si`` is composed of the following steps:

1. Create a MTF file set (Multiple TSV Files) describing your metabolic reactions and carbon transitions (.netw), experimental data (.miso) label input (.linp) and some non mandatory measurements options (.mflux, .mmet, .tvar, .cntsr, .opt). Let an example MTF set has a prefix ``mynetwork``. The syntax rules for reactions will be more or less obvious to someone working on metabolism biochemistry. So, to go quickly, you can inspire from example files ``test/mtf/e_coli.netw`` and co. distributed with the ``influx_si`` software (run ``influx_s --copy_test`` to bring them to your current directory). You can also consult the help message from ``txt2ftbl -h`` for ``--mtf`` option.

 .. note:: ``NA`` values (as "Non Available") are admitted as measurements values where appropriate. The difference with the situation where measurements are simply omitted is that NA measurements are simulated and are present in the vectors ``simulated unscaled labeling measurements`` and ``simulated scaled labeling measurements`` in the result kvh file.
 
 .. note:: In case of ``influx_i``, label kinetics can be provided in .miso file using non-empty ``Time`` column.
  Empty cells in ``Value`` column are equivalent to ``NA``.

2. Set your current directory to the directory of ``mynetwork.*`` and run ::

   $ influx_s.py --prefix mynetwork

  or ::

   $ influx_i.py --prefix mynetwork

  depending on stationary or instationary labeling context. We suppose here that directory of ``influx_si`` binaries is in the PATH variable.

  An ``influx_si`` run will produce the following files in the same directory that ``mynetwok.*``:
 
  ``mynetwork.ftbl``
    FTBL is a previously used format as a front-end format. It is still in use but behind the scenes. This file can be necessary as entry for ``ftbl2*`` utilities.
  ``mynetwork.log``
    contains the run-time output from various scripts, in particular,
    it contains a report on convergence history during the fitting process.
    It can be helpful for identifying potential problems, but if everything
    is going well, the user does not have to examine the content of this file;
  ``mynetwork.err``
   contains the warning and error messages.
   Normally, this file should be empty (0 byte size);
  ``mynetwork_res.kvh``
   contains all the results. `KVH format <http://serguei.sokol.free.fr/kvh-format/>`_ is a
   lightweight plain text format for hierarchically structured data. It can be seen in a text editor or in a spreadsheet software as its fields are tab separated. It can also be processed by user's custom software for post-processing, graphics output and alike. If ``influx_si`` is run on a series of starting points, there will be generated a common result
   file ``mynetwork_res.kvh`` which contains common information to all starting points
   but also a series of kvh files, one by starting point, e.g. ``mynetwork_res.V1.kvh``,
   ``mynetwork_res.V2.kvh`` and so on;
  ``mynetwork.pres.csv``
   contains a matrix of fitted parameters and final cost values. Each column
   corresponds to a particular starting point if run with ``--fseries`` and /or
   ``--iseries`` options. If ``influx_si`` was run without these options, the file
   will contain only one column corresponding to the starting point defined
   in the ``mynetwork.tvar`` file or to the random starting point.
  ``edge.netflux.mynetwok``, ``edge.xchflux.mynetwok``, ``node.log2pool.mynetwork``
   as the middle name of these files suggest, they can be used to map the corresponding
   values on the network graph in the `cytoscape <http://www.cytoscape.org>`_ software.
  
   .. note:: All these files are silently overwritten if already exist.
    So take care to copy your results elsewhere if you want to protect them
    from overwriting.
  
  custom files (e.g. ``mynetwork.pdf``)
    These files can be produced by user supplied scripts that are executed at the end of ``influx_si`` simulations. For example, we provide a script ``plot_ilab.R`` which can be used to plot label kinetics simulated by ``influx_i``. One or many of such custom scripts can be given in .opt file, field ``posttreat_R`` (cf. e_coli_i.opt for example)
  
  .. note:: It can be helpful to do some "dry runs" by executing ::

   $ influx_s.py --noopt --pref mynetwork
   
   before collecting actual measurement data to see if intended measurements will be sufficient to well define all fluxes, or at least the fluxes of interest. It is possible to do so because the measurement values in the .miso file have no impact on flux SD calculation when ``--noopt`` option is used. So it can be used any values, even NA at this moment. On the contrary, ``SD`` values set in .miso file, must be realistic. It is generally not a problem as they express measurements errors and are more or less known for a given measurement method.
   
   It is worthwhile to stress that a "dry run" is done for some presumed free flux values. If they reveal to be very different from actual flux values, it can happen that a network considered as well defined at moment of "dry run" turned into a badly defined network with actual measurement data and corresponding estimated fluxes. So it is important to do his best to guess the most realistic free fluxes for "dry runs" and log their values in .tvar file.

3. See warning and error messages in ``mynetwork.err`` if any. Correct what has to be corrected and retry p. 2

4. Extract and use the numerical results from the ``mynetwork_res.kvh`` file.

5. Optionally, visualize net fluxes (or exchange fluxes or logarithm of metabolite concentrations :math:`\log_2(M)`) in cytoscape using ``ftbl2xgmml`` to produce a .xgmml file and then mapping ``edge.netflux.mynetwok.attrs``, ``edge.xchflux.mynetwok.attrs`` or ``node.log2pool.mynetwork.attrs`` on graphical attributes like edge width, color etc. in cytoscape.
