.. _ftblevo:

.. highlight:: bash

=====================
FTBL format evolution
=====================

Even if FTBL is no more a front-end format for ``influx_si`` software (it was replace by MTF to this end), we consider useful to keep track of its evolution during ``influx_si`` development.


Introduction
------------
FTBL format was conceived by authors of ``13CFlux`` software in late 1990's (cf. https://www.13cflux.net/). At the beginning of 2000's, ``13CFlux`` became well spread in scientific community working on metabolism and isotope labeling. When we published the first version of ``influx_s`` in 2011, we adopted FTBL format to avoid cumbersome rewriting of networks and data already in use by the community. Second version of 13CFlux, published in 2012, abandoned FTBL format which was replaced by FluxML (XML) and was accompanied by a tool for automatic conversion of FTBL to FluxML.

On our side, we decided to continue to use FTBL by extending and evolving some of its features till its replacement by MTF (starting from v6.0). These extensions and evolution are presented hereafter for keeping tracks only. This chapter is not necessary for reading to successfully use ``influx_si`` software. Version number in titles indicates when described feature was first introduced to ``influx_si``.
 
``METABOLITE_POOLS`` and ``METAB_MEASUREMENTS`` (v2.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sections ``METABOLITE_POOLS`` and ``METAB_MEASUREMENTS`` concerning metabolite pools were added. These sections can be useful for stationary labeling when growth fluxes are modeled with :math:`\mu M` terms (cf. :ref:`Growth flux option <growthflux>`) or when some metabolites are confounded in measurements due to cell compartmentation of co-elution during HPLC step or whatever reason. These sections become mandatory for ``influx_i`` usage for instationary labeling as not only fluxes but also metabolite concentrations impact label propagation dynamics.

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Missing values marked as ``NA`` are admitted in measurement sections, in columns designated to values. In contrast, they are not admitted in columns designated to standard deviations. The main difference between a measurement just omitted and those marked as ``NA`` is that the latter will be simulated and reported in corresponding simulation sections of the result file.
This feature can be useful for preliminary simulations when there is no yet data available but user want to know e.g. if fluxes of interest will be well determined or not based on a supposed set of measurements. In this case, all presumed data can be set to ``NA`` (but not their SD).

Optimization control parameters (v5.3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Optimization method(s) can be tuned by control parameters that can be put in ``OPTIONS`` section. The format of those fields has changed. Before, the field names were looking like ``optctrl_maxit`` i.e. a prefix ``optctrl_`` followed by a parameter name, here ``maxit``. Starting from v5.3, they look like ``optctrl:nlsic:maxit`` i.e. a prefix ``optctrl`` followed by a method name (here ``nlsic``) and ended by parameter name, like ``maxit``, all 3 separated by colon ``:``. This new format allows tuning parameters for multiple optimization methods simultaneously. It became necessary, as starting from v5.3, several optimization methods can be used successively in one ``influx_si`` run. More about parameters can be found in the section :ref:`Optimization options <optopt>`.

