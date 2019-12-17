
.. _manual:

.. highlight:: bash

.. _MetExplore: http://metexplore.toulouse.inra.fr/

=============
User's manual
=============

Before diving in ``influx_si`` features let us present FTBL format evolution that was necessary to support ``influx_si`` innovations.

FTBL format evolution
---------------------
FTBL format was conceived by authors of ``13CFlux`` software in late 1990's (cf. https://www.13cflux.net/). At the beginning of 2000's, ``13CFlux`` became well spread in scientific community working on metabolism and isotope labeling. When we published the first version of ``influx_s`` in 2011, we adopted FTBL format to avoid cumbersome rewriting of networks and data already in use by the community. Second version of 13CFlux, published in 2012, abandoned FTBL format which was replaced by FluxML (XML) and was accompanied by a tool for automatic conversion of FTBL to FluxML.

On our side, we decided to continue to use FTBL by extending and evolving some of its features. These extensions and evolution are presented hereafter. Version number in titles indicates when presented feature was first introduced to ``influx_s(i)``.
 
``METABOLITE_POOLS`` and ``METAB_MEASUREMENTS`` (v2.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sections ``METABOLITE_POOLS`` and ``METAB_MEASUREMENTS`` concerning metabolite pools were added. These sections can be useful for stationary labeling when growth fluxes are modeled with :math:`\mu M` terms (cf. `Growth flux option`_) or when some metabolites are confounded in measurements due to cell compartmentation of co-elution during HPLC step or whatever reason. These sections become mandatory for ``influx_i`` usage for instationary labeling as not only fluxes but also metabolite concentrations impact label propagation dynamics.

``METABOLITE_POOLS`` is structured in two columns named ``META_NAME`` and ``META_SIZE`` and as ussual for FTBL indented and separated by tabulations, e.g. ::

		METABOLITE_POOLS
			META_NAME	META_SIZE
			AKG	-0.5
			...

.. note::
  
  The value ``-0.5`` is not aligned with its column name ``META_SIZE`` because by default, tab characters are expanded to 8 spaces. As ``META_NAME`` occupies 9 spaces, ``META_SIZE`` is just shifted to the next tab position. User has to use only one tab character to separate columns even if they don't look aligned on his screen.

For ``influx_i``, every internal metabolite (i.e. metabolites present in ``NETWORK`` section and not being input or output metabolites) and participating in carbon exchange must be referenced in this section. The value given in the column ``META_SIZE`` is a metabolite concentration. The unit used for these values must be in accordance with the units used for fluxes. For example, if metabolite concentrations are measured in mM/g then fluxes are supposed to be measured in mM/(g*[time_unit]). If the value is positive then corresponding metabolite is considered as having constant concentration which does not vary during fitting iterations. If the value is negative, then this metabolite concentration will be part of fitted variables and its absolute value is used as a starting value for these iterations.
A final fitted value will be expressed as a positive number.

For ``influx_s``, this section is optional and only few (not all) internal metabolites can be present in this section.

``METAB_MEASUREMENTS`` section regroups measurements of internal metabolite concentrations. Input and output metabolites may have concentrations varying during an experiment as they are consumed or produced. So they cannot appear in this section.  ``METAB_MEASUREMENTS`` section has 3 columns: ``META_NAME``, ``VALUE`` and ``DEVIATION``, e.g. ::

	METAB_MEASUREMENTS
		META_NAME	VALUE	DEVIATION
		Fru6P	0.43	0.01
		...

Column names are self explanatory.

In case of confounded measurements, confounded metabolites can be given as a sum, e.g. ::

	METAB_MEASUREMENTS
		META_NAME	VALUE	DEVIATION
		R5P_c+R5P_m	0.32	0.01
		...
  
In this case, the value ``0.32`` will be fitted by a sum of simulated metabolite concentrations.

Long reactions (v4.0)
~~~~~~~~~~~~~~~~~~~~~
Initially, FTBL admitted no more than 2 metabolites on each side of reactions put in ``NETWORK`` section. We had to overcome this limit to facilitate FTBL creation for studies including reactions much longer than that. Now, chemical reaction having more than two metabolites on any side can be split in several sub-reactions, each of which has no more then 2 metabolites on every side. It is important that all sub-reactions be put together one after another and that they  have the same name. Based on this name, ``influx_si`` will assemble all parts in one reaction. E.g. a reaction named ``Val_syn`` ::

  Val_syn: Pyr (abc) + Pyr (def) + Glu (ghijk) + NADPH -> Val (abcef) + CO2 (d) + AKG (ghijk)
  
can be translated into FTBL format as ::

 NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2
	Val_syn	Pyr	Pyr	Val	CO2
		#abc	#def	#abcef	#d
	Val_syn	Glu	NADPH	AKG	
		#ghijk	#	#ghijk	

If some reactions have the same name but not placed sequentially one after another, it will be signaled as an error.

Cofactors (v4.0)
~~~~~~~~~~~~~~~~~
Here, we call cofactors metabolites that does not participate in carbon transfer from one or several molecules to another. The main interest of entering cofactors in carbon transferring reactions is additional balance equations that we can put in stoechiometric system. Thus the number of free fluxes is diminished and fluxes are constrained to more realistic values, not violating cofactor balances.

To indicate that a metabolite is a cofactor, user can simply put an empty carbon string in the corresponding carbon transferring line. For example, a reaction ::

 v8: PEP (abc) -> Pyr (abc) + ATP
 
can be translated into FTBL as ::

 NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2
	v8	PEP		Pyr	ATP
		#abc		#abc	#

Note an empty carbon string ``#`` at the place corresponding to ``ATP``.
An important difference between cofactors and other metabolites that the former are allowed to have stoechiometric coefficients different from 1. These coefficients must be separated from cofactors by ``*`` sign, e.g. a reaction ::

  v41: Asp (abcd) + 2 ATP + NH3 -> Asn (abcd)

can be translated into FTBL as ::

 NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2
	v41	Asp	2*ATP	Asn	
		#abcd	#	#abcd	
	v41	NH3		
		#		

Note the presence of ``2*ATP`` term.

Same metabolite on both sides of reaction (v4.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some particular cases, it can be necessary to have a same metabolite on both sides of reaction. Let us illustrate this situation with the following example: ::

 v71: CO2.unlabeled (a) + CO2 (b) -> CO2 (a) + CO2.out (b)
 
Metabolite CO2 is present on both sides of reaction but its carbon atom is not the same. This is the main reason for introducing this feature, to allow tracer rearrangement. In FTBL, it gives ::

 NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2
	v71	CO2.unlabeled	CO2	CO2	CO2.out
		#a	#b	#a	#b


Section ``NOTRACER_NETWORK`` (v4.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to reactions with carbon rearrangements, it can be useful to add reactions with no carbon transfer. The most known reaction of such type is biomass composition but it can there be others, e.g. involving exclusively cofactors: ::

  v61: NADH + 0.5 O2 -> 2 ATP
  
This optional section is structured in 2 columns: ``FLUX_NAME`` and ``EQUATION``: ::

 NOTRACER_NETWORK
	FLUX_NAME	EQUATION
	v61	NADH+0.5*O2 = 2*ATP

You can see that the reaction is written in a manner very different form ``NETWORK`` section. Its sides are separated by ``=`` sign, metabolites are separated by ``+`` and they can have stoechiometric coefficients separated by ``*`` symbol. It is not visible in this example, but there can be as many metabolites as desired on each side of reaction. The limit "no more than 2 metabolites by side" proper to ``NETWORK`` section does not apply here.

Sub-sections ``EQUALITY/METAB`` and ``INEQUALITY/METAB`` (v2.11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the same manner as for fluxes, user can have to constrain variable metabolite concentrations. Constraints can be by equalities and inequalities. These subsections are organized in the same way as for fluxes. In ``EQUALITY/METAB`` there are 2 columns ``VALUE`` and ``FORMULA`` while in ``INEQUALITY/METAB`` there are 3 of them: ``VALUE``, ``COMP`` and ``FORMULA``. For example, ::

 EQUALITIES
	METAB
		VALUE	FORMULA
		0	R5P - 1.5*X5P
		...
 INEQUALITIES
	METAB
		VALUE	COMP	FORMULA
		0.001	<=	PEP
		10	>=	PEP
		...

``NA`` in measurements (v2.5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Missing values marked as ``NA`` are admitted in measurement sections, in columns designated to values. In contrast, they are not admitted in columns designated to standard deviations. The main difference between a measurement just omitted and those marked as ``NA`` is that the latter will be simulated and reported in corresponding simulation sections of the result file.
This feature can be useful for preliminary simulations when there is no yet data available but user want to know e.g. if fluxes of interest will be well determined or not based on a supposed set of measurements. In this case, all presumed data can be set to ``NA`` (but not their SD).

Convention evolution
~~~~~~~~~~~~~~~~~~~~~
Not only FTBL format evolved but also some conventions between its parts and content did so. Here is a complete list of them:
  - user must explicitly declare input-output fluxes as non reversible (set them as ``C`` with a value ``0`` in the section ``FLUX/XCH``) to make a distinction between input-output metabolites and "dead-end" metabolites (the latter are allowed since the version 2.0 and have net flux equal to 0 while exchange flux non zero).
  - starting from the version 2.8, new fluxes (i.e. absent in the ``NETWORK`` section) may appear in ``EQUALITY`` section. They can come, for example, from stoechiometry on cofactors involving non carbon carrying fluxes. These new fluxes have still to be declared in ``FLUX/{NET,XCH}`` sections (even if this feature is maintained in v4.0 its interest has diminished since cofactors can now be directly introduced in ``NETWORK`` and ``NOTRACER_NETWORK`` sections);
  - in LABEL_INPUT section following conventions apply since v3.2:
      * *"the rest is unlabeled"*: if many labeling forms are lacking in the file (including fully unlabeled metabolite) and the present forms does not sum up to 1, then the fully unlabeled form is considered as completing the set to 1;
      * *"guess the lacking one"*: if only one form is lacking in the file (no matter which one), then its fractions is considered as completing the present set to 1.
  - starting from v4.2, a particular comment tag ``//##`` is used to introduce a pathway name. The information on pathways can be useful for visualization on a partner web site MetExplore_ (cf. ``ftbl2metxml`` in _`Additional tools` section).

Basic influx_si usage
---------------------
``influx_si`` can be run without any option on most common cases. So its usage can be as simple as ::

 $ influx_s.py mynetwork
 
or ::

 $ influx_i.py mynetwork

we suppose here that a valid `FTBL <https://www.13cflux.net/>`_ file ``mynetwork.ftbl`` was created. Moreover, we supposed ``influx_s.py`` and ``influx_i.py`` is in the PATH variable.

In the rest of this manual, we'll use just ``influx_s.py`` as example if the example is valid for both stationary and instationary contexts. If some usage is valid exclusively for ``influx_i.py``, it will be duly signaled.


In a high throughput context, it can be useful to proceed many FTBL files in parallel. This can be done by giving all the FTBL names in a command line, e.g. ::

 $ influx_s.py mynetwork1 mynetwork2

and so on. All files are then proceeded in separate independent processes launched almost simultaneously by a bunch of size equal to the number of available or requested cores (if an option ``--np=NP`` is used). It is an operating system who is in charge to make a distribution of all these processes among all available CPUs and cores.

Sometimes, particular cases need usage of special options of ``influx_si``. The list of available options can be seen by running::

 $ influx_s.py --help

If used with options, ``influx_si`` can be run like ::

 $ influx_s.py [options] mynetwork

where ``[options]`` is an option list separated by a white character. Each option starts with a double dash ``--`` and can be followed by its argument if applicable. For example, to use BFGS optimization method instead of the default NLSIC algorithm, a user can run::

 $ influx_s.py --meth BFGS mynetwork

or ::

 $ influx_s.py --meth=BFGS mynetwork

The option names can be shortened till a non ambiguous interpretation is possible, e.g in the previous example, the option could be shortened as ``--m BFGS`` or ``--m=BFGS`` because there is no other option name starting by a letter ``m``. But an option ``--no`` could not be distinguished between ``--noopt`` and ``--noscale``. So at least ``--nos`` (for ``--noscale``) or ``--noo`` (for ``--noopt``) should be provided. There is only one option that does not admit a usage of an equal sign to provide an argument, it is ``--excl_outliers``. Use only a space character to provide an argument to this option when required.

Here after the available options with their full names are enumerated and detailed.

``influx_si`` command line options
----------------------------------
	--version        show program's version number and exit
	-h, --help       show the help message and exit
	--noopt          no optimization, just use free fluxes as is (after a projection on feasibility domain), to calculate
									 dependent fluxes, cumomers, stats and so on
	--noscale        no scaling factors to optimize => all scaling factors are assumed to be 1

									 This option can be useful if your measurements are already scaled to sum up to 1 which is often the case of MS data. Then, user saves some free parameters corresponding to scaling factors. This option can become mandatory if user wants to prevent scaling factors to be adjusted by optimization process.
	--meth=METH      method for optimization, one of nlsic|BFGS|Nelder-Mead.
									 Default: nlsic
	--fullsys        calculate all cumomer set (not just the reduced one
									 necessary to simulate measurements)

									 This option influences only post-optimization treatment. The fitting itself is still done with the reduced cumomer set or EMU variables if requested so. See the original paper on ``influx_s`` for more information on the reduced cumomer set.
	--emu            simulate labeling in EMU approach

									 This option should not produce a different result in parameter fitting. It is implemented and provided in a hope that on some network the results can be obtained in a shorter time
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

									 A byproduct of this option is that it can drastically reduce  cumomer system sizes. As it ensures that non reversible fluxes cannot change the sign, revers fluxes can be eliminated from pathways leading to observable cumomers. 
	--cinout=CINOUT  lower limit for input/output free and dependent fluxes.
									 Must be non negative. Default: 0
	--clowp=CLOWP    lower limit for free metabolite pools. Must be positive. Default 1.e-8
	--np=NP            When integer >= 1, it is a number of parallel threads (on
										 Unix) or subprocesses (on Windows) used in Monte-Carlo
										 (M-C) simulations or for multiple FTBL inputs. When NP is
										 a float number between 0 and 1, it gives a fraction of
										 available cores (rounded to closest integer) to be used.
										 Without this option or for NP=0, all available cores in a
										 given node are used for M-C simulations.
	--ln             Least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient

									 Jacobian can become rank deficient if provided data are not sufficient to resolve all free fluxes. It can be useful to determine fluxes that can still be resolved by the available measurements. If the Jacobian does not become rank deficient, this option has no influence on the found solution neither on the optimization process. But if the Jacobian does become rank deficient, a warning message is printed in the error file even if the optimization process could go to the end.

									 .. note:: Use this option with caution, in particular, when used in conjunction with Monte-Carlo simulations. As undetermined fluxes will be given some particular value, this value can be more or less stable from one Monte-Carlo simulation to another. This can create an illusion that a flux is well determined. See the linearized statistics in the result file to decide which fluxes are badly resolved.

									 A correct way to deal with badly defined metabolic network is to provide additional data that can help to resolve all the fluxes and/or to optimize input label, not just put ``--ln`` option and cross the fingers.

									 .. warning:: In this option, the notion of "least norm" is applied to *increments* during the optimization, not to the final solution. So undetermined fluxes could vary from one run to another if the optimization process is started from different points while well determined fluxes should keep stable values.
	--sln            Least norm of the solution of linearized problem (and not just of increments) is used when Jacobian is rank deficient
	--tikhreg        Approximate least norm solution is used for increments
									 during the non-linear iterations when Jacobian is rank
									 deficient
									 
									 To obtain an approximate solution a Tikhonov regularization is used when solving an LSI problem. Only one of the options ``--ln`` and ``--tikhreg`` can be activated in a given run.
	--lim            The same as --ln but with a function limSolve::lsei()
	--zc=ZC          Apply zero crossing strategy with non negative threshold
									 for net fluxes
									 
									 This option can accelerate convergence in situations when a net flux has to change its sign during the optimization iterations. Once such flux is identified, it is better to write the corresponding reaction in an opposite sens in the FTBL file or to give a starting value with a correct sign to avoid such zero crossing situation.
	--ffguess        Don't use free/dependent flux definitions from FTBL
									 file(s). Make an automatic guess.
									 
									 The fact that free fluxes are chosen automatically does not allow to specify a starting point for optimization iterations so a random starting point is used (drawn uniformly in [0; 1] interval). An option ``--seed`` can be useful to make the results reproducible.
	--fseries=FSERIES  File name with free parameter values for multiple
										 starting points. Default: '' (empty, i.e. only one
										 starting point from the FTBL file is used)
										 
										 The file must be formatted as plain text file with tab separator. There must be as many columns as starting points and at least as many rows as free parameters assigned in this file. A subset of free parameters can be used in this file. In this case, the rest of parameters take their unique starting values from the FTBL file. The first column must contain the names of free parameters used in this file. If there are extra rows whose names are not in the set of free parameter names, they are simply ignored. The first row must contain the names of starting points. These names can be just numbers from 1 to the number of starting points.
	--iseries=ISERIES  Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default '' (empty, i.e. all provided starting points are used)
										 
										 When used with conjunction with ``--fseries``, this option indicates the starting points to use from FSERIES file. But this option can also be used in conjunction with ``--irand`` to generate a required number of random starting points, e.g. ``influx_s.py --irand --iseries 1:10 mynetwork`` will generate and use 10 random starting points.
										 
										 For both ``--fseries`` and ``--iseries``, one result file is generated per starting point, e.g. ``mynetwork_res.V1.kvh``, ``mynetwork_res.V2.kvh`` and so on. If starting points comes from a ``--fseries`` then the suffixes ``V1``, ``V2``, ... are replaced by the column names from this file. In addition, a file ``mynetwork.pres.csv`` resuming all estimated parameters and final cost values is written.
	--seed=SEED        Integer (preferably a prime integer) used for
										 reproducible random number generating. It makes
										 reproducible random starting points (--irand) but also
										 Monte-Carlo simulations for sensitivity analysis.
										 Default: none, i.e. current system value is used, so
										 random drawing will be varying at each run.
	--excl_outliers    This option takes an optional argument, a p-value between
										 0 and 1 which is used to filter out measurement outliers.
										 The filtering is based on Z statistics calculated on
										 reduced residual distribution. Default: 0.01.

										 Excluded outliers (if any) and their residual values are reported in the ``mytework.log`` file. Non available (``NA``) measurements are considered as outliers for any p-value.
										 An optional p-value used here does not give a proportion of residuals that will be excluded from optimization process but rather a degree of beeing a valuable measurements. So, closer to zero is the p-value, the less data is filtered out. If in contary, you want to filter out more outliers than with the default p-value, use a value grater than the default value of 0.01, e.g.: ::

											influx_s.py --excl_outliers 0.02 mynetwork.ftbl

										 .. note::

											Don't use an equal sign "=" to give a p-value to this option. Here, only a white space can be used as a separator (see the example above).
	--nocalc          generate an R code but not execute it.
											
										This option can be useful for parallel execution of the generated R files via ``source()`` function in cluster environment
	--addnoise        Add centered gaussian noise to simulated measurements written to _res.kvh file. SD of this noise is taken from       FTBL file
	
										 This option can be helpful for generating synthetic FTBL files with realistic simulated measurements (cf. :ref:`How to make FTBL file with synthetic data?<howto>`).
	--TIMEIT          developer option

										Some portions of code are timed and the results is printed in the log-file. A curious user can use this option without any harm.
	--prof            developer option

										This option provides much more detailed profiling of the execution than ``--TIMEIT`` option. Only developers can be interested in using such information.

All command line options can be also provided in the FTBL file. A user can put them in the field ``commandArgs`` in the ``OPTIONS`` section. The corresponding portion of the FTBL file could look like

.. code-block:: none

 OPTIONS
	OPT_NAME	OPT_VALUE
	commandArgs	--meth BFGS --sens mc=100 --np 1

In such a way, a user can just drag-and-drop an FTBL file icon on the icon of the ``influx_s.py`` and the calculations will be done with the necessary options, assuming that the system was configured in appropriate way during the installation process.

If an option is provided both on the command line and in the FTBL file, it is the command line that has the priority. In such a way, a user is given an opportunity to overwrite any option at the run time. Nevertheless, there is no way to cancel a flag option (an option without argument) on a command line if it is already set in the FTBL file. For example, if ``--fullsys`` flag is set in the FTBL file, the full system information will be produced whatever command line options are.

Parallel experiments
--------------------
Staring from v4.0, ``influx_si`` offers possibility to treat labeling  data from parallel experiments. Parallel experiments for stationary labeling were described in the literature (e.g. cf. "Parallel labeling experiments and metabolic flux analysis: Past, present and future methodologies.", Crown SB, Antoniewicz MR., *Metab Eng.* 2013 Mar;16:21-32. doi: 10.1016/j.ymben.2012.11.010). But for instationary labeling, at the best of our knowledge, ``influx_si`` is the first software offering parallel experiments treatment.

The main interest of parallel experiments is increased precision of flux estimations. This comes at price of additional work for experiments and data gathering but the result is often worth the effort. As usual, before doing a real "wet" experiment, it can be useful to run few  "dry" simulations to see if planned experiments will deliver desired precision.

To deal with parallel experiments, a user have to prepare a series of FTBL files, one per experiment. One of them will be referred to as a main file. It has to provide the following sections common to all experiments: ``NETWORK``, ``FLUXES``, ``EQUALITIES`` (if any), ``INEQUALITIES`` (if any), ``FLUX_MEASUREMENTS`` (if any), ``METABOLITE_POOLS`` (if any), ``METAB_MEASUREMENTS`` (if any) and some entries in ``OPTIONS``. 

The secondary FTBL files as well as the main one are to provide experimental labeling data corresponding to each experiment. These data have to be presented in the following sections: ``LABEL_INPUT``, ``LABEL_MEASUREMENTS`` (if any), ``PEAK_MEASUREMENTS`` (if any), ``MASS_SPECTROMETRY`` (if any). In instationary context, text files with labeling kinetics have to be provided, one per experiment. Their names have to be placed in the field ``OPTION/file_labcin`` of a corresponding FTBL. Finally, the names of secondary FTBL files have to be put in the field ``OPTIONS/prl_exp`` of the main file as plain list separated by semicolon ``;`` and optionally by one or more spaces.

This file architecture ensures that a network topology, flux and metabolite values are common to all experiments while entry label and measurements on labeled metabolites are proper to each experiment.

Secondary FTBL files can also contain ``NETWORK`` and other sections found in the main file but are simply ignored at processing step. When FTBL files are ready, you can run ``influx_si`` on them by providing the name of main FTBL on the command line (and only it, don't list secondary files), e.g. in installation directory run: ::

 $ ./influx_s.py test/prl_exp/e_coli_glc1-6n

You can find an example of parallel experiment data in the directory ``test/prl_exp`` in files 
e_coli_glc1-6n.ftbl (main file), e_coli_glc2n.ftbl, e_coli_glc3n.ftbl, e_coli_glc4n.ftbl, e_coli_glc5n.ftbl, e_coli_glc6n.ftbl. These files correspond to stationary labeling experiments described in "Complete-MFA: Complementary parallel labeling experiments technique for metabolic flux analysis", Robert W. Leighty, Maciek R. Antoniewicz, *Metabolic Engineering* 20 (2013) 49–55 (with only difference that we use simulated and noised data instead of measured ones).

We also provide an example of simulated instationary parallel experiments in the files ``e_coli_GX_prl.ftbl`` (main file) and ``e_coli_GX_X.ftbl`` (secondary file) corresponding to simultaneous consumption of glucose and xylose. The network for this simulations was borrowed from "13C metabolic flux analysis of microbial and mammalian systems is enhanced with GC-MS measurements of glycogen and RNA labeling", Christopher P. Long, Jennifer Au, Jacqueline E. Gonzalez, Maciek R. Antoniewicz, Metabolic Engineering 38 (2016) 65–72. The experiment consisted in dynamic labeling by uniformly labeled glucose (main experiment)  and by uniformly labeled xylose (secondary one). Labeling kinetics MS data are given in ``e_coli_GX_MS.txt`` and ``e_coli_GX_X_MS.txt`` files respectively. To play with this example, you can run (still in installation directory): ::
 
 $ ./influx_i.py test/prl_exp/e_coli_GX_prl

The secondary files in all examples contain also the full information about the network, fluxes and so on, so they can be used as classical mono-experimental files to see how much the precision of flux estimation increased due to parallel experiment methodology.

Note that set of measured metabolite fragments as well as sampling time points for instationary labeling are not necessary the same for all parallel experiments. They do can differ.

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

	 btdesc=0.1
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
		return or not (default) the matrices with optimization steps and residual vectors during optimization. These matrices can then be found as part of ``optimization process information/history`` field in ``mynetwork_res.kvh`` file. Use it with caution, big size matrices can be generated requiring much of memory and disk space.

	 adaptbt=TRUE
		use (default) or not an adaptive backtracking algorithm.
		
	 monotone=FALSE
		should or not the cost decrease be monotone. If TRUE, then at first non decrease of the cost, the iterations are stopped with a warning message.

Names and default values for BFGS and Nelder-Mead algorithms can be found in the R help on ``optim()`` function.

Growth flux option
------------------
If present, this option makes ``influx_si`` take into account growth fluxes :math:`-\mu{}M` in the flux balance, where :math:`\mu` is a growth rate and :math:`M` is a concentration of an internal metabolite M by a unit of biomass. Only metabolites for which this concentration is provided in an FTBL section ``METABOLITE_POOLS``, contribute to flux balance with a flux :math:`-\mu{}M`.
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

Metabolite names used in this section must be identical to those used in the ``NETWORK`` section and others. Negative value is used as indicator of a variable metabolite pool. Such varying metabolites are part of fitted parameters. Absolute values from this section are used as their starting values in the optimization process.

One of valuable originality of ``influx_s``, it is a possibility to couple fluxomics and metabolomics in stationary experiments. It can be done because metabolite pools can influence labeling in two ways:

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

Like for other measurements, user has to provide a name, a value and a standard deviation for each entry in this section. Metabolites listed in this section must be defined in the ``NETWORK`` section and must have a negative value in the ``METABOLITE_POOLS`` section. Numerical values can be simple arithmetic expressions (as in the example above) which are evaluated during file parsing.

When a metabolite name is given as a sum of metabolites (e.g. ``Rub5P+Rib5P+Xul5P``) it is interpreted as a list of metabolites to be pooled. It is done proportionally to their concentrations. No numerical factor can appear in this sum. At least one of the metabolites from the list must have negative value in the ``METABOLITE_POOLS`` section. Otherwise, all metabolites from the list would be considered as having a fixed concentration and providing a measurement for such metabolites would be meaningless.

.. note:: There is no a specific option activating simulation of metabolite concentrations and taking them into account to the fitting process. Their simple presence in the ``METABOLITE_POOLS`` and ``METAB_MEASUREMENTS`` sections make concerned metabolites fittable parameters.

An example of an FTBL file having metabolite sections and involving growth fluxes can be found in ``test/e_coli_growth.ftbl``.

Post treatment option
---------------------

User can specify a name of one or several R scripts that will be automatically executed after non aborted ``influx_si`` run. This option can be useful, for example, for plain saving of calculation environment in a file for later exploring in an interactive R session or for plotting results in a pdf file and so on. A very basic example of such script is provided in the file ``test/save_all.R`` and its use can be found in the options of ``test/e_coli.ftbl`` file.

To activate this option, the script names must be provided in the ``OPTIONS`` section, in the field ``posttreat_R`` and separated by ``'; '``, e.g. ::

 OPTIONS
	OPT_NAME	OPT_VALUE
	posttreat_R	save_all.R; plot_something.pdf
	
The script name is interpreted as a relative path to the directory where the original FTBL file is located. After execution of ``save_all.R``, a file ``e_coli.RData`` is created. This particular example can be used to restore a calculation R environment by launching R and executing::

 > load("e_coli.RData")
 
After that, all variables defined in influx_si at the end of the calculations will be available in the current interactive session.
To be able to launch custom calculations on these variables, user has to do some preliminary actions. An example of such actions can be found in a file ``preamble.R`` which can be adapted for users's case.

To write his own scripts for post treatments or explore the calculated values in an interactive session, a user have to know some basics about existent variables where all the calculation results and auxiliary information are stored. Here are few of them:

dirw
	is a working directory (where the original FTBL file is)
dirx
	is an executable directory (where influx_s.py is)
baseshort
	is a short name of the input FTBL file (without the suffix ``.ftbl`` neither the directory part of the path)
param
	is the vector of the estimated parameters composed of free fluxes, scaling parameters (if any) and metabolite concentrations (if any)
jx_f
	is a environment regrouping calculated quantities. Here are some of its fields:
	
	fallnx
		a vector of all net and exchange fluxes (here, exchange fluxes are mapped on [0; 1[ interval)
	fwrv
		a vector of forward and reverse fluxes (reverse fluxes are "as is", i.e. not mapped)
	x
		is an internal state label vector
	simlab, simfmn and simpool
		are vectors of simulated measurements for label, net flux and metabolite pools respectively (fitting at the best of influx_s' capacity the provided measurements in the FTBL file)
	res
	 is the reduced residual vector, i.e. (simulated-measured)/SD
	ures
	 is the unreduced residual vector, i.e. (simulated-measured)
	jacobian
	 as its names indicates, is the Jacobian matrix (d res/d param)
	udr_dp
	 is the jacobian matrix for the unreduced residual vector (d ures/d param)

measurements
 is a list regrouping various measurements and their SD
nb_f
 is a list of various counts, like number of fluxes, parameters to fit, system sizes and so on
nm_list
 is a list of names for various vectors like fluxes, metabolites, label vectors, measurements, inequalities and so on
ui, ci
 are inequality matrix and right hand side respectively
 
A full list of all available variable and functions can be obtained in an R session by executing::

 > ls()
 
This list of more than 400 items is too long to be fully described here. We hope that few items succinctly described in this section will be sufficient for basic custom treatments.

An inspirations for your own custom treatments and/or plotting can be found in files ``plot_imass.R`` and ``plot_smeas.R`` that plot instationary and stationary data respectively in pdf files.

Exclusive ``influx_i`` options
------------------------------
There is only one exclusive option that can be given on a command line:

	--time_order=TIME_ORDER     Time order for ODE solving (1 (default), 2 or 1,2).
															Order 2 is more precise but more time consuming. The
															value '1,2' makes to start solving the ODE with the first
															order scheme then continues with the order 2.
															
															The scheme order can be important for the precision of flux and concentration estimations. The impact is not direct but can be very important. Please note that it can happen that order 1 fits the data with lower cost value function but it does not mean that the fluxes/concentrations are better estimated.

Other options occur as fields in the section ``OPTIONS`` of the FTBL file.

 ``file_labcin``
	 gives the name of the text file with label kinetics. If the file name starts with a "/", it is considered as 
	 
	 The values must be organized in a matrix where each row corresponds to a measured isotopomer/cumomer/mass-isotopologue while each column corresponds to a given time point. First column gives the names of labeled measured species and the first row contains time points.
	 
	 Matrix must be written one row per line and its entries (cells) must be separated by tabulations. Missing data can be signaled as ``NA`` or just an empty cell. Comments are allowed and must start with ``#`` sign. The rest of the line after ``#`` is simply ignored.
	 Empty lines are ignored. In such a way, comments can help to annotate the data and empty lines can help to format the file for better human readability.
	 All lines (a part from blank lines and comments) must have the same number of cells.
	 
	 The specie names must fit the names used in corresponding measurement sections of FTBL file. For example, a name ``m:Rib5P:1,2,3,4,5:0:693`` is composed of several fields separated by a column ``:``
	 
	 ``m``
		 indicates that data are of ``MASS_SPECTROMETRY`` type. Other possible values are ``l`` for ``LABEL_MEASUREMENTS`` and ``p`` for ``PEAK_MEASUREMENTS``
	 ``Rib5P``
		 metabolite name
	 ``1,2,3,4,5``
		 carbon numbers present in the measured fragment
	 ``0``
		 mass shift relative to fully unlabeled mass isotopologue: ``0`` corresponds to a fraction of unlabeled fragment, ``1`` to a fraction of fragments with only one labeled carbon atom and so on
	 ``693``
		 line number in FTBL file corresponding to this measurement. If previous fields are sufficient to unambiguously identify the measurement, this field can be omitted.
		 
	 Cf. ``test/e_coli_msne.txt`` (and corresponding ``test/e_coli_i.ftbl``) for more examples.
	 
	 The measurement precision (SD) is considered as constant during time and its values (one per measured specie) is given in the FTBL file, in the corresponding measurement section.
	 
	 All time points must be positive and put in increasing order. The time point 0 must be absent and is considered as labeling start. At that point all species are supposed to be fully unlabeled. This means also that all label measurements must be provided with a correction for natural 13C labeling. To prepare MS data with such correction, a software `IsoCor <https://metatoul.insa-toulouse.fr/metasys/software/isocor>`_ can help.
	 
	 There can be fictitious time points without any data in them. This feature can be used to increase the time resolution at some time intervals. The simulation of label propagation will be done and reported at these fictitious time points but the fitting will be obviously done only at time points having real data in them. For a regular time interval sub-division, it is more practical to use a parameter ``nsubdiv_dt`` (cf. hereafter) instead of fictitious time point in this file.
	 
	 If this field is empty or absent in the FTBL file then no fit can be done and a simple label simulation is calculated as if ``--noopt`` option were activated. Such simulation can be done only if a time grid is defined with the help of two other parameters: ``dt`` and ``tmax`` (cf. hereafter).
 ``nsubdiv_dt``
	 integer number of sub-intervals by which every time interval is divided to increase the precision of time resolution.
	 
	 It can happen that the value 1 (default) is sufficient for a satisfactory flux/concentration estimation. User can gradually increase this value (2, 3, ...) in successive ``influx_i`` runs to be sure that better time resolution does not impact parameter estimation. This property is called *grid convergence*. A grid convergence is necessary to overcome the result dependency on the choice of a numerical discretization scheme. A grid convergence can be considered as achieved when changes in estimated parameters provoked by a grid refinement are significantly lower than estimated confidence intervals for these parameters.
 ``dt``
	 a real positive number, defines a time step in a regular grid in absence of a file in ``file_labcin`` field.
	 If a file with label kinetics is well present then this parameter has no effect.
	 
	 A regular time grid for label simulations can be useful on preliminary stage when user only elaborates FTBL file and wants to see if label simulation are plausible. It can also help to produce simulated measurements (which can be extracted from the ``_res.kvh`` file) for further numerical experiments like studying convergence speed, parameter identifiability, noise impact and so on.
 ``tmax``
	 a real positive number, defines the end of a regular time grid if the field ``file_labcin`` is empty or absent. Parameters ``dt`` and ``tmax`` must be defined in such a way that there will be at least 2 time points greater then 0 in the time grid.
	 
	 If a file with label kinetics is well present then this parameter can be used to limit time grid on which simulations are done. If the value in ``tmax`` is greater then the maximal time value defined in the kinetics file then this parameter has no effect.
	 
.. note::
	It is very important that the values for time, flux and metabolite concentrations be expressed in concordant units. It would be meaningless to give time in minutes, fluxes in mM/h/g and concentrations in mM. This will lead to wrong results.
	
	For example, if the time is expressed in seconds and concentrations in mM/g then fluxes must be expressed in mM/s/g.
	
.. note::
	Option ``--noscale`` must be always activated for instationary calculations. So that for example, MS measurements must be always composed of fully measured fragments (i.e. with all isotopologues present) and normalized to sum up to 1.

Result file fields
------------------

Generally speaking, the names of the fields in the result KVH file are chosen to be self explanatory. So there is no so much to say about them. Here, we provide only some key fields and name conventions used in the result file.

At the beginning of the ``mynetwork_res.kvh`` file, some system information is provided. Here "system" should be taken in two sens: informatics and biological. The information is reported in the fields  ``influx`` and  ``system sizes``. These fields are followed by  ``starting point`` information regrouping ``starting free parameters``,  ``starting cost value``, ``flux system (Afl)`` and ``flux system (bfl)``. Name conventions used in these and other fields are following:

 net and exchange fluxes
	are prefixed by ``n.`` or ``x.`` respectively
 free, dependent, constrained and variable growth fluxes
	are prefixed by ``f.``, ``d.``, ``c.`` and ``g.`` respectively. So, a complete flux name could look like ``f.n.zwf`` which means `free net ZWF flux`.
	Growth fluxes which depend on constant metabolite concentrations can be found in constrained fluxes. Constant or variable growth fluxes are postfixed with ``_gr`` (as `growth`) string. For example, a flux ``g.n.Cit_gr`` corresponds to a net growth flux of Citrate metabolite. The growth fluxes are all set as non reversible, so all exchange fluxes like ``g.x.M_gr`` or ``c.x.M_gr`` are set to 0.
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

The field ``optimization process information`` is the key field presenting the results of an optimization process. The fitted parameters are in the subfield ``par``. Other subfields provide some additional information.

The final cost value is in the field ``final cost``.


The values of vectors derived from free fluxes like dependent fluxes, cumomers, MID and so on are in the corresponding fields whose names can be easily recognized.

Linear stats and Monte-Carlo statistics are presented in their respective fields. The latter field is present only if explicitly requested by user with ``--sens mc=MC`` option. In this kvh section, a term ``rsd`` means "relative standard deviation" (in literature, it is often encountered a synonym CV as Coefficient of Variation), it is calculated as SD/Mean and if expressed in percentage then the formula becomes 100%*SD/Mean.

The field ``jacobian dr_dp (without 1/sd_exp)`` report a Jacobian matrix which is defined as a matrix of partial derivatives :math:`\partial{r}/\partial{p}` where *r* is residual vector (Simulated--Measured) and *p* is a free parameter vector including free fluxes, scaling factors (if any) and free metabolite pools (if any). Note that in this definition the residual vector is not yet scaled by standard deviation of measurements. Sometimes, Jacobian is called *sensitivity matrix* in which case a special care should be brought to the sens of derivation. Often, by sensitivity matrix, we intend a matrix expressing how estimated fluxes are sensible to variations in the measurement data. Such definition corresponds to generalized inverse of Jacobian and it is reported in the field ``generalized inverse of jacobian dr_dp (without 1/sd_exp)``

Network values for Cytoscape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Several network values formatted for cytoscape are written by ``influx_si`` to their respective files. It can facilitate their visualizing and presentation in graphical mode. All these values can be mapped on various graphical attributes like edge width, node size or color scale of any of them. All these files are written at the end of calculations so if an error has interrupted this process, no such file will be produced. Take care to don't use an outdated copy of these files.

A file named ``edge.netflux.mynetwork.attrs`` can help to map net flux values on edges of a studied network. A file ``edge.xchflux.mynetwork.attrs`` do the same with exchange fluxes. And finally, ``node.log2pool.mynetwork.attrs`` provides logarithm (base 2) of pool concentrations. They can be mapped on some graphical attribute of network nodes.

See `Additional tools`_ section, `ftbl2xgmml: cytoscape view`_ paragraph to know how to produce files importable in Cytoscape from a given FTBL file. User's manual of Cytoscape has necessary information about using visual mapper for teaching how some values like net flux values can be mapped on graphical elements like edge width and so on.

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

	Error : Flux matrix is not square or singular: (56eq x 57unk)
	You have to change your choice of free fluxes in the 'mynetwork.ftbl' file.
	Candidate(s) for free flux(es):
	d.n.Xylupt_U

a message about badly structurally defined network could be similar to

.. code-block:: text

	Error : Provided measurements (isotopomers and fluxes) are not
		sufficient to resolve all free fluxes.
	Unsolvable fluxes may be:
		f.x.tk2, f.n.Xylupt_1, f.x.maldh, f.x.pfk, f.x.ta, f.x.tk1
	Jacobian dr_dff is dumped in dbg_dr_dff_singular.txt

a message about singular cumomer balance matrix could resemble to

.. code-block:: text

	lab_sim: Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3 or constrain some of the fluxes listed below to be non zero Zero rows in cumomer matrix A at weight 1:
	cit_c:16
	ac_c:2
	...
	Zero fluxes are:
	fwd.ACITL
	...


	
.. note:: In this error message, we report cumomers whose balance gave a zero row in the cumomer matrix (here ``cit_c:<N>`` cumomers, where <N> is an integer, its binary mask indicates the "1"s in the cumomer definition) as well as a list of fluxes having 0 value. This information could help a user to get insight about a flux whose zero value led to a singular matrix. A workaround for such situation could be setting in the FTBL file an inequality constraining a faulty flux to keep a small non zero value. A more radical workaround could be restricting some flux classes (input-output  fluxes with the option ``--cinout=CINOUT`` or even all non reversible ones with the option ``--clownr=CLOWNR``) to stay out of 0, e.g.:
 
 ``$ influx_s.py --clownr 0.0001 mynetwork``
 
 Adding such inequalities does not guaranty that cumomer matrix will become invertible but often it does help.
 It's up to user to check that an addition of such inequalities does not contradict biological sens of his network.

a message about badly statistically defined network could appear like

.. code-block:: text

 Inverse of covariance matrix is numerically singular.
 Statistically undefined parameter(s) seems to be:
 f.x.pyk
 For more complete list, see sd columns in '/linear stats'
 in the result file.

and so on.

A user should examine carefully any warning/error message and start to fix the problems by the first one in the list (if there are many) and not by the easiest or the most obvious to resolve. After fixing the first problem, rerun ``influx_si`` to see if other problems are still here. Sometimes, a problem can induce several others. So, fixing the first problem could eliminate some others. Repeat this process, till all the troubles are eliminated.

Problematic cases
-----------------

Obviously, everyone would like be able just run a flux estimation software and simply get results but unfortunately it does not work in this way every time.
In this section we review some problematic cases which can be encountered in practice.

Structurally non identifiable fluxes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can happen that collected data are not sufficient to resolve some fluxes in your network. Due to non linear nature of the problem, this situation can appear for some set of free flux values and disappear for others or be persistent for any free flux values. An error is reported to signal such situation, e.g.

.. code-block:: text

 lsi: Rank deficient matrix in least squares
 1 unsolvable variable(s):
 f.n.PPDK        7

and execution is stopped.

Several options are then available for a user facing such situation.

1. Collect more data to resolve lacking fluxes. As a rule of thumb, data must be collected on metabolites which are node of convergence of badly defined fluxes or on metabolites situated downhill of convergence point and preserving labeling pattern. Nature of collected data can be also important. Examples can be constructed where mass data are not sufficient to determine a flux but RMN data can do the job.
 
 Before using real data collection, you can make a "dry run" with ``--noopt`` option and with fictitious values for intended metabolite in the FTBL file to see if with these new data, the network becomes well resolved. If the error message disappear and SD values in the section ``linear stats`` are not very high then chances are that additionally collected data can help to resolve the fluxes.
 
2. Optimize input label. It can happen that you do collect data on a metabolite situated in convergence point for undefined fluxes but incoming fluxes are bringing the same labeling pattern which prevents flux(es) to be resolved. May be changing substrate label can help in this situation. For label optimization you can use a software called IsoDesign, distributed under OpenSource licence and available here http:://metatoul.insa-toulouse.fr/metasys/software/isodes/ (may be you have received ``influx_si`` as part of IsoDesign package, in which case you have it already).
 
 Naturally, this label optimization should be done before doing actual experiments. See IsoDesing tutorial for more details on how to prepare and make such optimization.
 
 If you don't want or don't have a possibility to use a software for label optimization or you think to have an insight on what should be changed in substrate labeling to better define the fluxes, you can still make a try with ``influx_s.py --noopt new_labeling.ftbl`` option to see if a new labeling will do the job (here ``new_labeling.ftbl`` is an example name for a FTBL file that you will prepare with a new ``LABEL_INPUT`` section.)

3. Use ``--ln`` option. It wont make you fluxes well defined, it will just continue calculation trying to resolve what can be solved and assigning some particular values (issued from so called *least norm* solution for rank deficient matrices) to undefined fluxes. You will still have a warning similar to

 .. code-block:: text

	 lsi_ln: Rank deficient matrix in least squares
	 1 free variable(s):
	 f.n.PPDK        7
	 Least L2-norm solution is provided.
 
 informing you that some flux(es) in the network is(are) still undefined. This option can be helpful if undefined fluxes are without particular interest for biological question in hand and their actual values can be safely ignored.

4. You can give an arbitrary fixed value to an undefined flux by declaring it as constrained in the FTBL file (letter ``C`` in the column ``FCD`` in the ``FLUXES`` section).

Badly defined fluxes
~~~~~~~~~~~~~~~~~~~~

Also known as *statistically undefined fluxes*, these fluxes have big or even huge SD values. The difference between these fluxes and structurally undefined fluxes is that the badly defined fluxes can become well defined if the noise is reduced or hypothetically eliminated while the latter will still be undetermined even in the absence of the noise. Despite this difference, all options presented in the previous section are applicable here (all but ``--ln`` which would be without effect here).

An additional measure can be taken which consist in experimental noise reduction. Generally, it can be done by using better protocols, better instruments or simply by increasing the measurement repetition number.

Once again, a use of ``--noopt`` with new hoped DEV values in the FTBL file can help to see if these new measurements with better noise characteristics will resolve or not the problem.

Slow convergence
~~~~~~~~~~~~~~~~

Slow optimization convergence can manifest by following warnings::

 nlsic: Maximal non linear iteration number is achieved

or/and ::

 nlsic: Maximal backtrack iteration number is achieved
 
Theoretically, user can increase the limit for those two numbers
(``optctrl_maxit`` and ``optctrl_btmaxit`` respectively in the ``OPTIONS`` section of FTBL file) but generally it is not a good idea. It can help only in very specific situations that we cannot analyze here as we estimate them low probable.
In all cases, a slow convergence is due to high non linearity of the solved problem. What can vary from one situation to another, it is the nature of this non linearity. Depending on this nature, several steps can be undertaken to accelerate optimization:

1. If a non linearity causing the slow convergence is due to the use of function absolute value :math:`|x|` in the calculation of forward and revers fluxes from net and exchange fluxes, then an option ``--zc=ZC`` (zero crossing) can be very efficient. This non linearity can become harmful when during optimization a net flux has to change its sign, in other words it has to cross zero.

 This option splits the convergence process in two parts. First, a minimum is searched for fluxes under additional constraints to keep the same sign during this step. Second, for fluxes that reached zero after the first step, a sign change is imposed and a second optimization is made with these new constraints.
 If ``--zc`` option is used with an argument 0 (``--zc=0`` or ``--zc 0``), it can happen that fluxes reaching zero produce a singular (non invertible) cumomer balance matrix. In this case, an execution is aborted with an error starting like
 
	.. code-block:: text
	 
		Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3 or constrain some of the fluxes listed below to be non zero
		...
	 
 To avoid such situation, an argument to ``--zc`` must be a small positive number, say ``--zc 0.001``. In this case, positive net fluxes are kept over 0.001 and negative fluxes are kept under -0.001 value. In this manner, an exact zero is avoided.
 
 Another way to avoid problem induced by using module function :math:`|x|` is to add inequality(-ies) imposing sens of reaction in ``INEQUALITIES/NET`` section, e.g. ::
	
	 0.0001	<=	mae
 
 Naturally, in this example, you have to be sure that the reaction catalyzed by malic enzyme (here ``mae``) must go in the sens written in your FTBL file.
 
 You can find potential candidates to impose sens of reaction by examining the flux values in ``mynetwork_res.kvh`` after a slow convergence and looking fluxes who's sign (positive or negative) looks suspicious to you. In our practice, we could observe a dramatic increase in convergence speed and stability just after imposing a sens of reaction to a "key" reaction. Obviously, such constraint must be in accordance with biological sens of a studied network and its biological condition.
 
2. A high non linearity can appear for some particular set of fluxes, especially when they take extreme values, e.g. when exchange fluxes are close to 1 or net fluxes take very high values of order 10² or even 10³ (supposing that the main entry flux is normalized to 1). In such a case, user can low this limits (options ``--cupx=CUPX`` and ``--cupn=CUPN`` respectively) or try to exclude outliers (``--excl_outliers P-VALUE``) as outliers can attract the solution in weird zone of fluxes. In this latter case, the first convergence will continue to be slow and will generate corresponding warnings but the second one (after a possible automatic elimination of outliers) can be much quicker.


Convergence aborted
~~~~~~~~~~~~~~~~~~~
This situation is signaled by an error message::

 nlsic: LSI returned not descending direction

This problem can occur for badly defined network which are very sensitive to truncation errors. The effect of such errors can become comparable to the effect of the increment step during optimization. It means that we cannot decrease the norm of residual vector under the values resulting from rounding errors.
If it happens for relatively small increments then the results of convergence are still exploitable. If not, there is no so many actions that user could undertake except to make his system better defined as described in previous sections.

.. note:: By default, we use a very small value for increment norm as stopping criterion (:math:`10^{-5}`). It can be considered as very drastic criterion and can be relaxed to :math:`10^{-3}` or :math:`10^{-2}` depending on required precision for a problem in hand (to do that, use an option ``optctrl_errx`` in the section ``OPTIONS`` of FTBL file). 

Additional tools
----------------

Tools described in this section are not strictly necessary for running ``influx_si`` and calculating the fluxes. But in some cases, they can facilitate the task of tracking and solving potential problems in FTBL preparation and usage.

Most of the utilities produce an output written on standard output or in a file who's name is derived from the input file name. This latter situation is signaled with a phrase "The output redirection is optional" and in the usage examples the output redirection is taken in square brackets ``[> output.txt]`` which obviously should be omitted if an actual redirection is required. Such behavior is particularly useful for drag-and-drop usage.

txt2ftbl: conversion of txt format to FTBL format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An easily readable/writable text format can be used to create *de novo* an FTBL file. Reactions in this text format can look like:

.. code-block:: text
 
	v48: Asp (abcd) + Pyr (efg) + Glu (hijkl) + SucCoA (mnop) + ATP + 2 NADPH ->
	 LL-DAP (0.5 abcdgfe + 0.5 efgdcba) + AKG (hijkl) + Suc (0.5 mnop + 0.5 ponm)

This long reaction illustrates several format features:

 - ``v48`` is the reaction name. It is optional. If reaction names (and theirs separators ``:`` signs) are omitted, reactions will be just numbered. The numbering restarts after each comment block (a comment starts with a ``#`` sign). This is done to give an opportunity to organize reactions in pathways. In such a case, a comment is considered as stating a new pathway which is also numbered. Thus an automatic reaction name can look like ``r2.3`` where ``2`` is pathway number and ``3`` is a reaction number in this pathway.
 
  Note that in this example, the reaction is split in two lines only for convenience of presentation. In a text file, a reaction must be written on only one line. No line breaks are admitted in reactions and no more than one reaction can be written on a given line.
  
 - ``Asp (abcd)`` is a metabolite name ``Asp`` followed by its optional carbon id string between parentheses ``(abcd)``. All carbon id must be a unique letter on each side of the reaction and if present on one side of reaction, must also be present on the other one. Thus carbon atom balance is preserved. In case of symmetric molecule, the carbon scrambling can be indicated as ``(0.5 mnop + 0.5 ponm)`` as e.g. for succinate ``Suc`` in the example above. Numeric coefficients of carbon forms (here ``0.5``) can be omitted as all forms are considered as equally probable and automatically normalized to sum up to 1. So a completely equivalent form could be ``(mnop + ponm)``.
 - ``+`` sign separates metabolites on each side of reaction
 - ``->`` separates two sides of reaction and indicates that this reaction is irreversible, i.e. its exchange flux is zero. It does not precludes about the sens of reaction. Here we consider that a reaction can be irreversible and have a negative net flux. If in addition, you wish to indicate that a reaction must operate only from left to right, i.e. to have a positive net flux, then use ``->>`` sign. To indicate a reversible reaction use ``<->`` and a reversible reaction with imposed positive net flux use ``<->>``.
 - ``ATP`` is an example of a cofactor, it does not have a carbon id string. It participates in mass balance but not in carbon balance equations.
 - ``2 NADPH`` is an example of a cofactor with a stoechiometric coefficient different from 1. Coefficients different from 1 are not allowed for metabolites participating in carbon exchanges in a given reaction. But if a reaction has no carbon exchanges, then all metabolites are allowed to have a coefficient different from 1 like for example in biomass reactions.
 - ``#`` starts a comment that will be put in FTBL as is, except the first ``#`` hash sign that will be replaced by ``//`` (FTBL comment tag)
 - ``###`` triple hash sign is used to introduce a pathway name. Respectively, ``//##`` will do the same in FTBL. Pathway name can be useful for ``ftbl2metxml.py`` script which prepare xml and txt files for visualization on a partner site MetExplore_.
 
An example of a full featured metabolite network can be found in ``test/prl_exp/e_coli_anto.txt``.

To convert it to FTBL file, you can run (in installation directory):

 .. code-block:: text

	$ ./txt2ftbl.py test/prl_exp/e_coli_anto.txt [> test/prl_exp/e_coli_anto.ftbl]
 
Note that output redirection ``> ...`` is optional. In absence of such redirection, the output file name is guessed from input file by replacing ``.txt`` with ``.ftbl`` extension.
Thus obtained FTBL file must be completed with several kinds of information like label input, label measurements and so on to be fully functional and suitable for ``influx_si``.

ftbl2xgmml: cytoscape view
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once a valid FTBL file is generated, a user can visualize a graph representing his metabolic network in `Cytoscape <http://www.cytoscape.org>`_ program. To produce necessary graph files, user can run::

 $ ftbl2xgmml.py mynetwork[.ftbl] [> mynetwotk.xgmml]

or drag and drop ``mynetwork.ftbl`` icon on ``ftbl2xgmml.py`` icon.

The output redirection is optional.

This will produce a file in the XGMML format ``mynetwork.xgmml`` in the directory of ``mynetwork.ftbl``:

Once a generated file ``mynetwork.ftbl`` is imported in cytoscape, a user can use one of automatic cytoscape layouts or edit node's disposition in the graph by hand.
For those who use `CySBML <http://apps.cytoscape.org/apps/cysbml>`_ plugin, a saving of a particular layout in a file can be practical for later applying it to a new network.

Graphical conventions used in the generated XGMML are the following:

* metabolite are presented as rounded square nodes;
* simple (one to one) reaction are represented by simple edges;
* condensing and/or splitting reactions are represented by edges converging and/or diverging from additional almost invisible node having a label with the reaction name;
* all nodes and edges have tool tips, i.e. when a pointer is put over, their name (metabolite or reaction) appears in a tiny pop-up window;
* non reversible reaction are represented by a single solid line, have an arrow on the target end (i.e. produced metabolite) and nothing on the source end (i.e. consumed metabolite);
* reversible reactions are represented by a double parallel line and have a solid circle on the source end;
* color code for arrows:

	* green for free net flux;
	* blue for dependent net flux;
	* black for constrained net flux;

* color code for solid circles:

	* green for free exchange flux;
	* blue for dependent exchange flux;
	* black for constrained exchange flux.

ftbl2netan: FTBL parsing
~~~~~~~~~~~~~~~~~~~~~~~~

To see how an FTBL file is parsed and what the parsing module "understands" in the network, a following command can be run::

 $ ftbl2netan.py mynetwork[.ftbl] [> mynetwork.netan]

The output redirection is optional.

A user can examine ``mynetwork.netan`` in a plain text editor (not like Word) or in spreadsheet software. It has an hierarchical structure, the fields are separated by tabulations and the field values are Python objects converted to strings.

ftbl2cumoAb: human readable equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, it can be helpful to examine visually the equations used by ``influx_si``. These equations can be produced in human readable form by running::

 $ ftbl2cumoAb.py -r mynetwork[.ftbl] [> mynetwork.sys]

or::

 $ ftbl2cumoAb.py --emu mynetwork[.ftbl] [> mynetwork.sys]
 
The output redirection is optional.

The result file ``mynetwork.sys`` will contain systems of stoichiometric and cumomer balance equations as well as a symbolic inversion of stoichiometric matrix, i.e. dependent fluxes are represented as linear combination of free and constrained fluxes and an optional constant value. In the examples above, the option ``-r`` stands for "reduced cumomer set" and ``--emu`` stands for "generate EMU framework equations". In this latter case, only isotopologues of mass+0 in each EMU are reported in ``mynetwork.sys`` file. For other mass weights, equations does not change and the right hand side term could get longer for condensation reactions but involves the same EMUs as in mass+0 weight.

If a full cumomer set has to be examined, just omit all options. Keep in mind that on real-world networks this can produce more than thousand equations by cumomer weight which could hardly be qualified as *human* readable form. So use it with caution.

For the sake of brevity, cumomer names are encoded in decimal integer form. For example, a cumomer ``Metab#xx1x`` will be referred as ``Metab:2`` because a binary number ``0010`` corresponds to a decimal number ``2``. The binary mask ``0010`` is obtained from the cumomer mask ``xx1x`` by a plain replacement of every ``x`` by ``0``.

For a given cumomer weight, the equations are sorted alphabetically.

expa2ftbl: non carbon carrying fluxes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some reactions of carbon metabolism require cofactor usage like ATP/ADP and some others. A mass balance on cofactors can produce additional useful constraints on the stoechiometric system. Since the version 2.8, such mass balance equation on non carbon carrying metabolites can be put in ``EQUATION`` section of FTBL file. A utility ``expa2ftbl.R`` can be helpful for this purpose if a user has already a full set of reactions in `expa <http://gcrg.ucsd.edu/Downloads/ExtremePathwayAnalysis>`_ format.
To extract additional equation from an expa file, ``expa2ftbl.R`` can be used as::

 $ R --vanilla --slave --args file.expa < expa2ftbl.R > file.ftbl_eq

Then an information for the generated ``file.ftbl_eq`` has to be manually copy/pasted to a corresponding FTBL file.

Note that ``expa2ftbl.R`` uses a Unix command ``grep`` and another utility described here above ``ftbl2netan.py``.

res2ftbl_meas: simulated data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During preparation of a study, one of questions that biologist can ask is "Will the intended collected data be sufficient for flux resolution in a given network?"
Some clue can be obtained by making "dry runs" of ``influx_si`` with ``--noopt`` (i.e. no optimization) option. User can prepare an FTBL file with a given network and supposed data to be collected. At first, the measurement values can be replaced by NAs while the SD values for measurements must be given in realistic manner. After running::

 $ influx_s.py --noopt mynetwork

a utility ``res2ftbl_meas.py`` can be practical for preparing FTBL files with obtained simulated measurements::

 $ res2ftbl_meas.py res2ftbl_meas.py mynetwork_res[.kvh] > mynetwork.ftbl_meas

(here ``.kvh`` suffix is optional). The information from the generated file ``mynetwork.ftbl_meas`` has to be manually copy/pasted into corresponding FTBL file.
Getting an ftbl file with real values instead of NAs in measurement sections gives an opportunity to explore optimization behavior near a simulated point like convergence speed and/or convergence stability to cite few of them.

ffres2ftbl: import free fluxes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This utility imports free flux values and metabolite concentrations (if any) from a result file _res.kvh and inject them into an FTBL file. Usage::

 $ ffres2ftbl.sh mynetwork_res.kvh [base.ftbl] > new.ftbl

If an optional argument ``base.ftbl`` is omitted, then the free flux values are injected into an FTBL file corresponding to the _res.kvh file (here ``mynetwork.ftbl``). This script can be used on a Unix (e.g. Linux, MacOS) or on a cygwin (unix tools on Windows) platform. It makes use of another utility written in python ``ff2ftbl.py``

ftbl2kvh: check ftbl parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This utility simply parses a ftbl file and write what was "understood" in a kvh file. No network analysis occurs here unlike in ``ftbl2netan`` utility. Usage::

 $ ftbl2kvh.py mynetwork[.ftbl] [> mynetwork.kvh]

The output redirection is optional.

ftbl2metxml: prepare MetExplore_ visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Convert an FTBL file to an xml file suitable for visualization on MetExplore_ site. If a result kvh file ``mynetwork_res.kvh`` is present, it will be parsed to extract flux values corresponding to the last ``influx_si`` run and put them in ``mynetwok_net.txt``, ``mynetwork_fwd.txt`` and ``mynetwork_rev.txt``. As their names indicate, they will contain net, forward and revers flux values respectively.

IsoDesign: optimizing input label
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of means to increase a flux resolution can be an optimization of input label composition. A utility ``IsoDesing`` solving this problem was developed by Pierre Millard. It is not part of ``influx_si`` distribution and can be downloaded at http://metatoul.insa-toulouse.fr/metasys/software/isodes/. In a nutshell, it works by scanning all possible input label compositions with a defined step, running ``influx_si`` on each of them then collecting the SD information on all fluxes for all label compositions and finally selecting an input label optimal in some sens (according to a criterion chosen by a user).

.. _Cytoscape: http://www.cytoscape.org
