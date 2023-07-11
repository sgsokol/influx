
.. _quick:

===========
Quick Start
===========

A basic work-flow with ``influx_si`` is composed of the following steps:

1. Create a MTF file set (Multiple TSV Files) describing your metabolic reactions and carbon transitions (.netw), experimental data (.miso) label input (.linp), some non mandatory measurements (.mflux, .mmet) and optional files providing initial values for estimated parameters, constraints and so on (.tvar, .cntsr, .opt). Let an example MTF set have a prefix ``mynetwork``. The syntax rules for reactions will be more or less obvious to someone working on metabolism biochemistry. So, to go quickly, you can inspire from example files ``test/mtf/e_coli.netw`` and co. distributed with the ``influx_si`` software (run ``influx_s --copy_test`` to bring them to your current directory). You can also consult the help message from ``txt2ftbl -h`` for ``--mtf`` option.

 .. note:: ``NA`` values (as "Non Available") are admitted as measurements values where appropriate. The difference with the situation where measurements are simply omitted from a file is that NA measurements are simulated and are present in the result files ``*.sim`` while absent measurements are not.
 
 .. note:: In case of ``influx_i``, label kinetics can be provided in .miso file using non-empty ``Time`` column.
  Empty cells in ``Value`` column are equivalent to ``NA``.

2. Set your current directory to the directory of ``mynetwork.*`` and run ::

   $ influx_s.py --prefix mynetwork

  or ::

   $ influx_i.py --prefix mynetwork

  depending on stationary or instationary labeling context. We suppose here that directory of ``influx_si`` binaries is in the PATH variable.

  An ``influx_si`` run will produce result files in ``mynetwok_res`` directory. The detailed description of these files can be found in the next section. However, general idea is that simulated measurements are written in the files similar to MTF format with ``.sim`` (as "simulated") suffix appended, e.g. ``mynetwork.miso.sim``.
    
  .. note:: It can be helpful to do some "dry runs" by executing ::

   $ influx_s.py --noopt --pref mynetwork
   
   before collecting actual measurement data to see if intended measurements will be sufficient to well define all fluxes, or at least the fluxes of interest. It is possible to do so because the measurement values in the .miso file have no impact on flux SD calculation when ``--noopt`` option is used. So it can be used any values, even NA at this moment. On the contrary, ``SD`` values set in .miso file, must be realistic. It is generally not a problem as they express measurements errors and are more or less known for a given measurement method.
   
   It is worthwhile to stress that a "dry run" is done for some presumed free flux values. If they reveal to be very different from actual flux values, it can happen that a network considered as well defined at moment of "dry run" turned into a badly defined network with actual measurement data and corresponding estimated fluxes. So it is important to do his best to guess the most realistic free fluxes for "dry runs" and log their values in .tvar file.

3. See warning and error messages in ``mynetwork.err`` if any. Correct what has to be corrected and retry p. 2

4. Extract and use the numerical results from the ``mynetwork_res/*.sim`` files.

5. Optionally, visualize net fluxes (or exchange fluxes or logarithm of metabolite concentrations :math:`\log_2(M)`) in cytoscape using ``ftbl2xgmml`` to produce a .xgmml file and then mapping files from ``mynetwork_res/tmp`` (``edge.netflux.mynetwok.attrs``, ``edge.xchflux.mynetwok.attrs`` or ``node.log2pool.mynetwork.attrs``) to graphical attributes like edge width, color etc. in cytoscape.
