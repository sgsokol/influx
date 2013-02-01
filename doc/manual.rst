
.. _manual:

=============
User's manual
=============

``influx_s`` is designed to run without any option on most common cases. So its usage can be as simple as::

 $ influx_s.py mynetwork

we suppose here that a valid `FTBL <https://www.13cflux.net/>`_ file ``mynetwork.ftbl`` was created

.. note::
 A documentation on FTBL syntax rules can be found in its original place, i.e. in the documentation on 13CFlux software freely available at https://www.13cflux.net/

moreover we supposed ``influx_s.py`` is in the PATH variable.


Sometimes, particular cases need usage of special options of ``influx_s``. The list of available options can be seen by running::

 $ influx_s.py --help

If used with options, ``influx_s`` can be run like::

 $ influx_s.py [options] mynetwork

where ``[options]`` is an option list separated by a white character. Each option starts with a double dash ``--`` and can be followed by its argument if applicable. For example, to use BFGS optimization method instead of the default NLSIC algorithm, a user can run::

 $ influx_s.py --meth BFGS mynetwork

or::

 $ influx_s.py --meth=BFGS mynetwork

The option names can be shortened till a non ambiguous interpretation is possible, e.g in the previous example, the option could be shortened as ``--m BFGS`` or ``--m=BFGS`` because there is no other option name starting by a letter ``m``. But an option ``--no`` could not be distinguished between ``--noopt`` and ``--noscale``. So at least ``--nos`` (for ``--noscale``) or ``--noo`` (for ``--noopt``) should be provided.

Here after the available options with their full names are enumerated and detailed.

Command line options
--------------------
  --version        show program's version number and exit
  -h, --help       show the help message and exit
  --noopt          no optimization, just use free fluxes as is, to calculate
                   dependent fluxes, cumomers, stats and so on
  --noscale        no scaling factors to optimize => all scaling factors are assumed to be 1

                   This option can be useful if your measurements are already scaled to sum up to 1 which is often the case of MS data. Then, the user saves some free parameters corresponding to scaling factors. This option can become mandatory if the user wants to prevent scaling factors to be adjusted by optimization process.
  --meth=METH      method for optimization, one of nlsic|BFGS|Nelder-Mead.
                   Default: nlsic
  --fullsys        calculate all cumomer set (not just the reduced one
                   necessary to simulate measurements)

                   This option influences only post-optimization treatment. The fitting itself is still done with the reduced cumomer set. See the original paper on ``influx_s`` for more information on the reduced cumomer set.
  --emu            simulate labeling in EMU approach

                   This option should not produce a different result in parameter fitting. It is implemented and provided in the hope that on some network the results can be obtained in a shorter time
  --irand          ignore initial approximation for free parameters (free fluxes and metabolite concentrations) from the FTBL file or from a dedicated file (cf --fseries and --iseries
                   option) and use random values drawn uniformly from [0,1]
                   
                   It is recommended to use this option in conjunction with "--zc 0" option.
  --sens=SENS      sensitivity method: SENS can be 'mc[=N]', mc stands for
                   Monte-Carlo. N is the number of Monte-Carlo simulations.
                   Default for N: 10

                   The sensitivity information (i.e. the influence of the noise in the data on the estimated parameter variation) based on linearized statistics is always provided. So the user has to use this option only if he wants to compare this linearized information to the Monte-Carlo simulations. Note that the default value 10 for the number of simulations is far from to be sufficient to get reliable statistical estimations. This default option allows only to quickly check that this option is working as expected.
  --cupx=CUPX      upper limit for reverse fluxes. Must be in interval [0, 1]. Default: 0.999
  --cupn=CUPN      upper limit for net fluxes. Default: 1.e3
  --cupp=CUPP      upper limit for metabolite pool. Default: 1.e5
  --clownr=CLOWNR  lower limit for not reversible free and dependent fluxes.
                   Zero value (default) means no lower limit
  --cinout=CINOUT  lower limit for input/output free and dependent fluxes.
                   Must be non negative. Default: 0
  --clowp=CLOWP    lower limit for free metabolite pools. Must be positive. Default 1.e-8
  --np=NP          Number of parallel process used in Monte-Carlo simulations
                   Without this option or for NP=0 all available cores in a
                   given node are used

                   At the time of this writing, a third-part R module ``multicore``, on which we are based for parallel execution, is not stable enough on Windows platform. So a Windows user should use this option with argument 1, e.g.
                   ``influx_s.py --sens mc=100 --np 1 mynetwork``
  --ln             Approximate least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient

                   Jacobian can become rank deficient if provided data are not sufficient to resolve all free fluxes. It can be useful to determine fluxes that can still be resolved by the available measurements. If the Jacobian does not become rank deficient, this option has no influence on the found solution neither on the optimization process. But if the Jacobian does become rank deficient, a warning message is printed in the error file even if the optimization process could go to the end.

                   .. note:: Use this option with caution, in particular, when used in conjunction with Monte-Carlo simulations. As undetermined fluxes will be given some particular value, this value can be more or less stable from one Monte-Carlo simulation to another. This can create an illusion that a flux is well determined. See the linearized statistics in the result file to decide which fluxes are badly resolved.

                   A correct way to deal with badly defined metabolic network is to provide additional data that can help to resolve all the fluxes, not just put ``--ln`` option and cross the fingers.

                   .. warning:: The notion of "least norm" is applied to increments during the optimization, not to the final solution. So undetermined fluxes could vary from one run to another if the optimization process is started from different points while well determined fluxes should keep stable values.
  --zc=ZC          Apply zero crossing strategy with non negative threshold
                   for net fluxes
                   
                   This option can accelerate convergence in situations when a net flux has to change its sign during the optimization iterations. Once such flux is identified, it is better to write the corresponding reaction in an opposite sens in the FTBL file or to give a starting value with a correct sign to avoid such zero crossing situation.
  --fseries=FSERIES  File name with free parameter values for multiple
                     starting points. Default: '' (empty, i.e. only one
                     starting point from the FTBL file is used)
                     
                     The file must be formated as plain text file with tab separator. There must be as many columns as starting points and at least as many rows as free parameters assigned in this file. A subset of free parameters can be used in this file. In this case, the rest of parameters take their unique starting values from the FTBL file. The first column must contain the names of free parameters used in this file. If there are extra rows whose names are not in the set of free parameter names, they are simply ignored. The first row must contain the names of starting points. These names can be just numbers from 1 to the number of starting points.
  --iseries=ISERIES  Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default '' (empty, i.e. all provided starting points are used)
                     
                     When used with conjunction with ``--fseries``, this option indicates the starting points to use from FSERIES file. But this option can also be used in conjunction with ``--irand`` to generate a required number of random starting points, e.g. ``influx_s.py --irand --iseries 1:10 mynetwork`` will generate and use 10 random starting points.
                     
                     For both ``--fseries`` and ``--iseries``, one result file is generated per starting point, e.g. ``mynetwork_res.V1.kvh``, ``mynetwork_res.V2.kvh`` and so on. If starting points comes from a ``--fseries`` then the suffixes ``V1``, ``V2``, ... are replaced by the column names from this file. In addition, a file ``mynetwork.pres.csv`` resuming all estimated parameters and final cost values is written.
  --seed=SEED        Integer (preferably a prime integer) used for reproducible random number generating. It makes reproducible random starting points (``--irand``) but also Monte-Carlo simulations for sensitivity analysis (``--sens mc=N``) if executed in sequential way (``--np=1``). Default: current system value, i.e. the random drawing will be varying at each run.
  --DEBUG          developer option

                   Produce a lot of run-time information in the log-file and many additional files. This also can slow down the program in a drastic way. Don't use this option unless your know what your are doing.
  --TIMEIT         developer option

                   Some portions of code are timed and the results is printed in the log-file. A curious user can use this option without any harm.
  --prof           developer option

                   This option provides much more detailed profiling of the execution than ``--TIMEIT`` option. Only developers can be interested in using such information.

All command line options can be also provided in the FTBL file. A user can put them in the field ``commandArgs`` in the ``OPTIONS`` section. The corresponding portion of the FTBL file could look like

.. code-block:: none

 OPTIONS
	OPT_NAME	OPT_VALUE
	commandArgs	--meth BFGS --sens mc=100 --np 1

In such a way, a user can just drag-and-drop an FTBL file icon on the icon of the ``influx_s.py`` and the calculations will be done with the necessary options, assuming that the system was configured in appropriate way during the installation process.

If an option is provided both on the command line and in the FTBL file, it is the command line that has the priority. In such a way, a user is given an opportunity to overwrite any option at the run time. Nevertheless, there is no way to cancel a flag option (an option without argument) on a command line if it is already set in the FTBL file. For example, if ``--fullsys`` flag is set in the FTBL file, the full system information will be produced whatever command line options are.

Optimization options
--------------------
These options can help to tune the convergence process of the NLSIC (or any other chosen algorithm). They can be given only in an FTBL file, in the section OPTIONS. These options are prefixed with ``optctrl_`` which is followed by a particular option name. For example, ``optctrl_errx`` corresponds to the stopping criterion hereafter and the corresponding FTBL portion could look like

.. code-block:: none

 OPTIONS
	OPT_NAME	OPT_VALUE
	optctrl_errx	1.e-3

All possible options and their default values for NLSIC algorithm follow:

   errx=1.e-5
    stopping criterion. When the L2 norm of the increment vector of free parameters is below this value, the iterations are stopped.

   maxit=50
    maximal number for non-linear iterations.

   btstart=1.
    backtracking starting coefficient

   btfrac=0.25
    backtracking fraction parameter. It corresponds to the alpha parameter in the paper on ``influx_s``

   btdesc=0.75
    backtracking descending parameter. It corresponds to the beta parameter in the paper on ``influx_s``

   btmaxit=15
    maximal number of backtracking iterations

   trace=1
    report (=1) or not (=0) minimal convergence information

   rcond=1.e10
   condition number over which a matrix is considered as rank deficient

   ci=list(p=0.95, report=F)
    confidence interval reporting. This option is own to ``nlsic()`` function. It has no impact on the reporting of linear stats information in the result kvh file after the post-optimization treatment. This latter is always done.

   history=FALSE
    return or not (default) the matrices with optimization steps and residual vectors during optimization. These matrices can then be found as part of ``optimization process informations/history`` field in ``mynetwork_res.kvh`` file. Use it with caution, big size matrices can be generated requiring much of memory and disk space.

   adaptbt=TRUE
    use (default) or not an adaptive backtracking algorithm.

Names and default values for BFGS and Nelder-Mead algorithms can be found in the R help on ``optim()`` function.

Growth flux option
------------------
If present, this option makes ``influx_s`` take into account growth fluxes :math:`-\mu{}M` in the flux balance, where :math:`\mu` is a growth rate and :math:`M` is a concentration of an internal metabolite M by a unit of biomass. Only metabolites for which this concentration is provided in an FTBL section ``METABOLITE_POOLS``, contribute to flux balance with a flux :math:`-\mu{}M`.
This flux can be varying or constant during optimization process depending on whether the metabolite M is part of free parameters to fit or not. Usually, taking into account of this kind of flux does not influence very much on the estimated flux values. So, this option is provided to allow a user to be sure that it is true in his own case.

The option is activated by a field ``include_growth_flux`` in the ``OPTIONS`` section:

.. code-block:: none

 OPTIONS
	OPT_NAME	OPT_VALUE
	include_growth_flux	1

Value 0 cancels the contribution of the growth fluxes to the general flux balance.

Another necessary option is ``mu`` giving the value of `µ`:

.. code-block:: none

 OPTIONS
	OPT_NAME	OPT_VALUE
	mu	0.12

Finally, the metabolite concentrations by a unit of biomass are reported in a section ``METABOLITE_POOLS`` as:

.. code-block:: none

 METABOLITE_POOLS
	META_NAME	META_SIZE
	Fum	2.47158569399681
	Suc	-15.8893144279264
	Mal	-6.47828321758155
	...	...

Metabolite names used in this section must be identical to those used in the NETWORK section and others. Negative value is used as indicator of variable metabolite pools. Such varying metabolites are part of fitted parameters. Absolute values from this section are used as starting values in optimization process.

One of valuable originality of ``influx_s``, it is a possibility given to
the user to couple fluxomics and metabolomics in stationary experiments. It can be done because metabolite pools can influence labeling in two ways:
 * through metabolite pooling (due to compartmentation and/or coelution during chromatography)
 * through growth fluxes.

This last influence is often of low intensity compared to metabolite transformation fluxes. In literature, it is often neglected.

.. note:: ``METABOLITE_POOLS`` section was not present in the original FTBL format. It is added `ad hoc` and it is possible that its presence makes fail other software using such FTBL.

Another section that was added "ad hoc" to FTBL file is ``METAB_MEASUREMENTS``:

.. code-block:: none

 METAB_MEASUREMENTS
	META_NAME	VALUE	DEVIATION
	Suc	15.8893144279264*1.e-3/10.7	1.e-2
	Mal	6.47828321758155*1.e-3/10.7	1.e-2
	Rub5P+Rib5P+Xul5P	1.66034545348219*1.e-3/10.7	1.e-2

Like for other measurements, user has to provide a name, a value and a standard deviation for each entry in this section. Metabolites listed in this section must be defined in the ``NETWORK`` section and must have a negative value in the ``METABOLITE_POOLS`` section. Numerical values can be simple arithmetic expressions which are evaluated during file parsing.

When a metabolite name is given as a sum of metabolites (e.g. ``Rub5P+Rib5P+Xul5P``) it is interpreted as a list of metabolites to be pooled. It is done proportionally to their concentrations. No numerical factor can appear in this sum. At least one of the metabolites from the list must have negative value in the ``METABOLITE_POOLS`` section. Otherwise, all metabolites from the list would be considered as having a fixed concentration and providing a measurement for such metabolites would be meaningless.

An example of an FTBL file having metabolite sections and involving growth fluxes can be found in ``test/e_coli_growth.ftbl``.

Result file fields
------------------

Generally speaking, the names of the fields in the result KVH file are chosen to be self explanatory. So there is no so much to say about them. Here, we provide only some key fields and name conventions used in the result file.

At the beginning of the ``mynetwork_res.kvh`` file some system information is provided. Here "system" should be taken in two sens: informatics and biological. The informations are reported in the fields  ``influx`` and  ``system sizes``. These fields are followed by  ``starting point`` information regrouping ``starting free parameters``,  ``starting MID vector`` (MID stands for Mass Isotopomer Distribution),  ``starting cumomer vector``, forward-revers fluxes, net-exchange fluxes, starting residuals and some other subfields. Name conventions used in these and other fields are following:

 net and exchange fluxes
  are prefixed by ``n.`` or ``x.`` respectively
 free, dependent and constrained fluxes
  are prefixed by ``f.``, ``d.`` and ``c.`` respectively. So, a complete flux name could look like ``f.n.zwf`` which means `free net ZWF flux`.
 scaling factors names
  are formed according to a pattern similar to ``label;Ala;1`` which corresponds to the first group of measurements on Alanine molecule in labeling experiments. Other possible types of experiments are ``peak`` and ``mass``.
 MID vector names
  are looking like ``METAB+N`` where ``METAB`` is metabolite name and ``N`` goes from 0 to the number of carbon atoms in the considered molecule.
 cumomer names
  follow classical convention ``METAB#pattern_of_x_and_1``, e.g. ``Ala#x1x``
 forward and reverse fluxes
   are prefixed by ``fwd.`` and ``rev.`` respectively, e.g. ``fwd.zwf`` or ``rev.zwf``
 measurement names
   have several fields separated by a colon ``:``. For example, ``l:Asp:#xx1x:694`` deciphers like:

     * ``l`` stands for `labeling` experiment (others possibilities are ``p`` for `peak`, ``m`` for `mass` and ``pm`` for `metabolite pool`)
     * ``Asp`` is a metabolite name
     * ``#xx1x`` is a measurement identification
     * ``694`` is a line number in the FTBL file corresponding to this measurement.

The field ``optimization process informations`` is the key field presenting the results of an optimization process. The fitted parameters are in the subfield ``par``. Other subfields provide some additional informations.

The final cost value is in the field ``final cost``.


The values of vectors derived from free fluxes like dependent fluxes, cumomers, MID and so on are in the corresponding fields whose names can be easily recognized.

Linear stats and Monte-Carlo statistics are presented in their respective fields. The latter field is present only if explicitly requested by the user with ``--sens mc=MC`` option.

Network values for Cytoscape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Several network values formatted for cytoscape are written by ``influx_s`` to their respective files. It can facilitate their visualizing and presentation in graphical mode. All these values can be mapped on various graphical attributes like edge width, node size or color scale of any of them. All these files are written at the end of calculations so if an error has interrupted this process, no such file will be produced. Take care to don't use an outdated copy of these files.

A file named ``edge.netflux.mynetwork`` can help to map net flux values on edges of a studied network. A file ``edge.xchflux.mynetwork`` do the same with exchange fluxes. And finally, ``node.log2pool.mynetwork`` provides logarithm (base 2) of pool concentrations. They can be mapped on some graphical attribute of network nodes.

See `Additional tools`_ section, `Cytoscape view`_ paragraph to know how to produce files importable in Cytoscape from a given FTBL file. User's manual of Cytoscape has necessary information about using visual mapper for teaching how some values like net flux values can be mapped on graphical elements like edge width and so on.

Warning and error messages
--------------------------
The warning and error messages are logged in the ``.err`` suffixed file. For example, after running::

 $ influx_s mynetwok

the warnings and errors will be written in the ``mynetwork.err`` file.
This kind of messages are important for user not only to be aware that during calculations something went wrong but also to understand what exactly went wrong and to have an insight on how to fix it.

Problems can appear in all stages of a software run:

* parsing FTBL files
* R code writing
* R code execution

  * vector-matrix initialization
  * optimization
  * post-optimization treatment

Most of the error messages are automatically generated by underlying languages Python and R. These messages can appear somewhat cryptic for a user unfamiliar with these languages. But the most important error messages are edited to be as explicit as possible. For example, a message telling that free fluxes are badly chosen could look like::

  Error : Flux matrix is not square: (56eq x 57unk)
  You have to change your choice of free fluxes in the 'mynetwork.ftbl' file.
  Candidate(s) for free flux(es):
  d.n.Xylupt_U
  Execution stopped

a message about badly structurally defined network could be similar to::

  Error : Provided measures (isotopomers and fluxes) are not
    sufficient to resolve all free fluxes.
  Unsolvable fluxes may be:
    f.x.tk2, f.n.Xylupt_1, f.x.maldh, f.x.pfk, f.x.ta, f.x.tk1
  Jacobian dr_dff is dumped in dbg_dr_dff_singular.txt
  Execution stopped

a message about singular cumomer balance matrix could resemble to::

  Error in solve(A, b) : 
    cs_lu(A) failed: near-singular A (or out of memory)
  Error in trisparse_solv(lAb$A, lAb$b, iw, method = "sparse") : 
    Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3
  or constrain some of the fluxes listed below to be non zero
  Zero rows in cumomer matrix A at weight 1:
  PHB:4
  PHB:1
  PHB:2
  PHB:8
  Zero fluxes are:
  fwd.AACOAR_1
  fwd.ACOAAT
  ...
  Calls: opt_wrapper -> nlsic -> r -> param2fl_x -> trisparse_solv
  Execution stopped
  
.. note:: In this error message we report cumomers whose balance gave a zero row in the cumomer matrix (here ``PHB:<N>`` cumomers, where <N> is an integer, its binary mask indicates the "1"s in the cumomer definition) as well as a list of fluxes having 0 value. This information could help a user to get insight about a flux whose zero value led to a singular matrix. A workaround for such situation could be setting in the FTBL file an inequality constraining a faulty flux to keep a small non zero value. A more radical workaround could be restricting some flux classes (input-output  fluxes with the option ``--cinout=CINOUT`` or even all non reversible ones with the option ``--clownr=CLOWNR``) to stay out of 0, e.g.:
 
 ``$ influx_s.py --clownr 0.0001 mynetwork``
 
 Adding such inequalities does not guaranty that cumomer matrix will become invertible but often it does help.
 It's up to user to check that an addition of such inequalities does not contradict biological sens of his network.

a message about badly statistically defined network could appear like::

 Inverse of covariance matrix is numerically singular.
 Statistically undefined parameter(s) seems to be:
 f.x.pyk
 For more complete list, see sd columns in '/linear stats'
 in the result file.

and so on.

A user should examine carefully any warning/error message and start to fix the problems by the first one in the list (if there are many) and not by the easiest or the most obvious to resolve. After fixing the first problem, rerun ``influx_s`` to see if other problems are still here. Sometimes, a problem can induce several others. So, correcting the first problem could eliminate some others. Repeat this process, till all the troubles are eliminated.

Additional tools
----------------

Tools described in this section are not strictly necessary for running ``influx_s`` and calculate the fluxes. But in some cases, they can facilitate the task of tracking and solving potential problems in FTBL preparation and usage.

Cytoscape view
~~~~~~~~~~~~~~

Once a valid FTBL file is generated, a user can visualize a graph representing his metabolic network in `Cytoscape <http://www.cytoscape.org>`_ program. To produce necessary graph files, user can run::

 $ ftbl2rsif.py mynetwork

or drag and drop ``mynetwork.ftbl`` icon on ``ftbl2rsif.py`` icon.

This will produce a series of files in the directory of ``mynetwork.ftbl``:

 .. describe:: mynetwork.sif

   this file has to be imported in Cytoscape (File > Import > Network (Multiple File Types)...)

 .. describe:: edge.targetArrowShape.mynetwork

 .. describe:: edge.targetArrowColor.mynetwork

 .. describe:: edge.sourceArrowShape.mynetwork

 .. describe:: edge.sourceArrowColor.mynetwork

 .. describe:: edge.label.mynetwork

   these files define graphical attributes of edges and should be imported via ``File > Import > Edge Attributes ...``
 .. describe:: node.shape.mynetwork

 .. describe:: node.fillColor.mynetwork

   these files define node visual attributes and should be imported via ``File > Import > Node Attributes ...``

Once all import finished, a user can use one of automatic cytoscape layouts or edit node's disposition in the graph by hand.

FTBL parsing
~~~~~~~~~~~~

To see how an FTBL file is parsed and what the parsing module "understands" in a given FTBL, a following command can be run::

 $ ftbl2netan.py mynetwork > mynetwork_netan.kvh

The end part of the command ``> mynetwork_netan.kvh`` means that the standard output (typically a console display) will be redirected to a file named ``mynetwork_netan.kvh``. A user can examine this file which has an hierarchical structure and where the values are Python objects converted to strings.

Human readable equations
~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, it can be helpful to examine visually the equations used by ``influx_s``. These equations can be produced in human readable form by running::

 $ ftbl2cumoAb.py -r mynetwork > mynetwork.sys

The result file ``mynetwork.sys`` will contain systems of stoichiometric and cumomer balance equations as well as a symbolic inversion of stoichiometric matrix, i.e. dependent fluxes are represented as linear combination of free and constrained fluxes and an optional constant value. In the example above, the option ``-r`` stands for "reduced cumomer set". If a full cumomer set has to be examined, just omit ``-r`` option. Keep in mind that on real-world networks this can produce more than thousand equations by cumomer weight which could hardly be qualified as *human* readable form. So use it with caution.

For the sake of brevity, cumomer names are encoded in decimal integer form. For example, a cumomer ``Metab#xx1x`` will be referred as ``Metab:2`` because a binary number ``0010`` corresponds to a decimal number ``2``. The binary mask ``0010`` is obtained from the cumomer mask ``xx1x`` by a plain replacement of every ``x`` by ``0`` .

For a given cumomer weight, the equations are sorted alphabetically.

An option ``--emu`` will generate symbolic equations for EMU framework instead of cumomer ones. Only isotopologues of mass+0 in each EMU are reported in this file. For other mass weights, equations does not change and the right hand side term could gets longer for condensation reactions but involves the same EMUs as in mass+0 weight.

.. _Cytoscape: http://www.cytoscape.org
