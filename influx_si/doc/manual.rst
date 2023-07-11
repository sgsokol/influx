
.. _manual:

.. highlight:: bash

.. _MetExplore: https://metexplore.toulouse.inra.fr/
.. _vkvh: https://vkvh.readthedocs.io
.. _Cytoscape: https://www.cytoscape.org

=============
User's manual
=============

.. _mtf:

MTF format
~~~~~~~~~~

Let start by describing file formats used by ``influx_si``.
MTF (Multiple TSV files) format was introduced in v6.0 of ``influx_si`` and was developed in collaboration and on demand from `MetaToul-FluxoMet platform <https://www.toulouse-biotechnology-institute.fr/en/plateformes-plateaux/metatoul/>`_ and `MetaSys team <https://www.toulouse-biotechnology-institute.fr/en/poles/equipe-metasys/>`_ of Toulouse Biotechnology Institute. Since v7.0, MTF is also used for output files.

Originally, MTF introduction came in replacement to FTBL format and served several purposes:

 - to simplify network and other information formatting;
 - to allow multiplexing of constant and variable information for large experience series. E.g., a same network can be tested under several biological conditions or vice versa a given data set can be confronted to different network models to see which one can better fit it;
 - to facilitate automatic assembling of network/data/options coming from tiers workflows such as MS or NMR data treatment.
 
 Even if FTBL format was replaced as front-end format in ``influx_si``, it is still a valid
 format for calculation and continues to be used behind the scene in accompanying utilities.
 
MTF format is composed of a series of plain text (tab separated when necessary) files, each having a particular role and designated by their suffixes:

 .netw
   describes biochemical reactions and label transitions.
 .linp
   describes label inputs forms and fractions
 .miso
   describes stationary and instationary label measurements, their type, value, and standard deviation
 .mflux
   describes stationary net flux measurements (e.g., substrate consumption or product secretion)
 .mmet
   describes stationary concentration measurements
 .cnstr
   describes constraints, equality or inequality type, on fluxes and concentrations
 .tvar
   describes variable type such as Free, Constrained or Dependent with possible starting values for Free variables
 .opt
   describes options such as optimization parameters, post-treatment R scripts to execute and others

These suffixes will be used throughout this manual, as they facilitate ``influx_si`` usage. But they are not mandatory. For example, a classical .tsv or .txt suffixes can be used instead. In this case, a file meaning can be indicated in ``--mtf`` option of ``influx_si`` of ``txt2ftbl``, e.g. ``--mtf netw=ecoli.txt,miso=gcms.tsv,linp=glucoseU.8.tsv``

Only the first 3 are mandatory as input files, the rest of MTF set is optional.

Now, we pass in more detailed review each of the file content.

.netw
-----

This is a central file to the rest of the work -- network file .netw
It represents a list of biochemical reactions with label transitions. Here is an example of such reaction: ::

 edd: Gnt6P (ABCDEF) -> Pyr (ABC) + GA3P (DEF)

Let detail its elements (all reaction elements are case-sensitive):
 ``edd``
   the unique reaction name
 ``:``
  name-reaction separator
 ``Gnt6P``, ``Pyr``, ``GA3P``
   reactant names. They must not contain ``+``, ``(``, ``)``. Other common separators such as ``:``, ``;`` and special characters such as ``{``, ``|``, ...  are to avoid. However, other than Latin alphabet letters are welcome in UTF-8 encoding, e.g. ``α``, ``β`` etc.
 ``(ABCDEF)``, ``(ABC)``, ``(DEF)``
   are labels of reactants. Here, their meaning is that first 3 carbons of ``Gnt6P`` goes to ``Pyr`` in the same order and the rest goes to ``GA3P`` also in the same order. Label atom numbering is left free for user's choice but once chosen, it must remain consistent between reactions. It is advised to follow some common conventions. E.g., the most oxidized carbon atom (if it is a carbon which is used for labeling) should have number 1 (here letter ``A`` in ``Gnt6P`` and ``Pyr``, and letter ``D`` in ``GA3P``). Its neighbor has number 2 etc. ``influx_si`` does not impose to use ¹³C as labeling atom. User is free to use any other atom or even a combination of them, e.g., ¹³C and ¹⁵N. A reactant can be without label atoms. In this case, it will participate in mass balance but not in label balance. Such situation can be useful e.g., for co-factors, e.g., ATP/ADP, NADP/NADPH etc. if they are not synthesized in the modeled network. The left opening parenthesis of the label pattern must be separated by a space from a precedent reactant name. The label identifiers can be composed of Latin or foreign letters as well as digits, each symbol representing one label atom. If a label symbol is present on one reaction side, it must be present as well on the other side. Each symbol must be present only once, except for so called "scrambling" label molecules. For example, Formate is a symmetric molecule from carbon point of view. So, its labeling can be given as ``(ABCD+DCBA)``. Thus, for example, a letter ``A`` appearing on the other side of reaction can come from both ends of Formate.
 ``->``
   reaction left and right side separator. Here, the reaction is indicated as non-reversible, i.e., exchange flux is imposed to be 0. However, the net flux can be either positive or negative. For a positive net flux, the reactants on the left are consumed and those on the right are produced. If user wishes to impose the sens of reaction, a sign ``->>`` can be used. In this case, the net flux is imposed to be non-negative (i.e. ≥ 0). For reversible reactions, a sign ``<->`` can be used. A sign ``<->>`` is meaningful. It indicates that a reaction is reversible, but the net flux must be non-negative. A non-zero exchange flux for a positive net flux can be responsible for a backward label propagation.
  
   Input/output reactions must have exchange flux set to 0. So they can be edited either with ``->`` or ``->>``. This is necessary to distinguish them from reactions having so-called "dead-end" reactants. They occur only on one side of reaction(s) (either left or right, like input/output reactants do) but have an exchange flux different from 0. The dead-ends are rarely desirable and most often result from network topology errors or simply leaving ``<->`` where ``->`` is meant. However, there are some special situations when they can be used on purpose in ``influx_i``, for example, for modeling label dilution from stock species.
 ``+``
   reactant separator. The surrounding spaces are mandatory.
 
 In the above reaction example, all stoichiometric coefficients are 1. If it is not the case in some reaction, they can be given as plain number (integer or with decimal point) preceding a reactant name and separated by an optional ``*`` sign (it can be replaced by a space), e.g.: ::
 
  v47:	Ser (abc) + AcCoA (de) + 3.0*ATP () + 4.0*NADPH () + SO4 () -> Cys (abc) + Ac (de)
 
 Note that coefficients different from 1 can only be used with non labeled reactants.
 
 A comment can be introduced by ``#``. The line content starting from this character to the end of the line is simply ignored with one exception: a triple hash sign ``###`` at the line beginning is used to introduce a pathway name. Pathway name can be useful for ``ftbl2metxml.py`` script which prepare xml and txt files for visualization on a partner site MetExplore_.

 
.linp
-----

Label input can be indicated in .linp file. Starting from this file type, the rest of the files are in TSV (tab separated values) format. I.e. they are plain text files where data are organized in tables, one table row per file line and where columns are separated with the tabulation character. The first non commented row contains column names. The comments start with ``#`` sign, they are simply ignored till the end of the row where they occur. The left- and right-trailing white spaces are stripped.

The .linp file can contain the following column names:
  ``Id``
    Not used by ``influx_si`` but can be useful for user's information tracking system.
  ``Comment``
    Not used by ``influx_si`` but left for user's convenience.
  ``Specie``
    Chemical specie name such as used in reactions in .netw file, e.g. ``Glucose``
  ``Isotopomer``
    Particular isotopomer form used as label entry and composed only of "0"s and "1"s, e.g. ``111111`` for uniformly labeled Glucose, ``100000`` for Glucose labeled in first carbon.
  ``Value``
    For step-wise labeling experiments, a number between 0 and 1 indicating a fraction of the given labeled form in the mix. For experiments with labeling varying in time, a R expression describing time dependent function for the given label form.
    
    Normally, all labeling form for a given specie must sum up to 1. If they don't, the following conventions apply:
    
    * *"the rest is unlabeled"*: if many labeling forms are lacking in the file (including fully unlabeled specie) then the fully unlabeled form (e.g. ``000000`` for Glucose) is considered as completing the set to 1;
      
    * *"guess the lacking one"*: if only one form is lacking in the file (no matter which one), then its fractions is considered as completing the present set to 1.
    
    Here is a complete example: ::
    
	Id	Comment	Specie	Isotopomer	Value
			Gluc_U	111111		1
			Gluc_U	000000		0.
			Gluc_1	100000		1.
			Gluc_1	000000		0.

    which can be shortened, due to above conventions, to::
    
	Id	Comment	Specie	Isotopomer	Value
			Gluc_U	111111		1
			Gluc_1	100000		1.
    

.miso
-----

The .miso contains and describes labeled measurements.
This file contains the following column names:

  ``Id``
    Not used by ``influx_si`` but can be useful for user's information tracking system.
  ``Comment``
    Not used by ``influx_si`` but left for user's convenience.
  ``Specie``
    Specie names such as used in reactions in .netw file, e.g. ``Glucose``
  ``Fragment``
    Integer sequences or intervals describing label fragment used in a given measurement, e.g. ``1-3`` or ``1,2,3`` or ``2-5,7,9-11``. Empty field means that the entire molecule is measured.
  ``Dataset``
    Any character sequence identifying measurement method. It can be useful for distinguishing measurements on the same combination specie/fragment. Examples: ``MS-1``, ``HSQC``. A given dataset can have its own scaling factor if they are in use.
  ``Isospecies``
    For MS measurements, a ``M0``, ``M1`` etc. isotopologue identification. For NMR label measurements, a combination of binary cumomers involved in measurements. They are separated by "+" sign. Each binary cumomer can be composed of "0", "1" and "x" symbols, e.g. ``01x+10x``. This notation is universal enough to describe any NMR (or even MS) method. However, for methods focused on a particular specie fragment, it can be more practical to use notation "label transferring" like ``2->``, ``2->1``, ``2->3`` and ``2->1,3``. Here, the second atom is labeled and in a given measurement method it interacts:
    
      ``2->``
          with no other atom (given a singlet peak);
      ``2->1``
          with labeled atom 1 (giving a peak doublet 1);
      ``2->3``
          then labeled atom 3 (giving peak doublet 2);
      ``2->1,3``
          and finally with both labeled atoms 1 and 3 (doublet of doublets).
      
  ``Value``
    Measured value in floating point notation. Can be empty or NA, meaning "non-available". There is a difference between a measurement absent in file and a measurement with NA value. The former is simply ignored, while the latter is simulated and reported is simulated measurements.
  ``SD``
    Standard deviation value in floating point notation. Cannot be empty, neither NA. We recall that SD is characterizing a given measurement technique, not its particular realization. So it is perfectly possible to have only one measurements if SD of the given technique was already estimated from previous experiments. If the ``Value`` contains an average of :math:`n` measurements, than a standard SD should be reduced by a factor :math:`\sqrt{n}`
  ``Time``
    For instationary labeling, the time point to which a given measurement corresponds. For stationary labeling, must be empty. 

A multi-line example is the following: ::

	Id	Comment	Specie	Fragment	Dataset	Isospecies	Value		SD	Time
			GA3P			LAB-10	1xx		0.03304418	0.002	
			GA3P			LAB-11	x1x		0.01260362	0.002	
			GA3P			LAB-12	xx1		0.1207158	0.002	
			PEP	1,2		PEAK-1	1->		0.66		0.005	
			PEP	1,2		PEAK-1	1->2		0.01		0.005	
			PEP	1,2,3		PEAK-2	2->		1.26		0.005	
			PEP	1,2,3		PEAK-2	2->1		0.03960991	0.005	
			PEP	1,2,3		PEAK-2	2->3		0.01183004	0.005	
			PEP	1,2,3		PEAK-2	2->1,3		0.0007686513	0.005	
			PEP	1,2,3		MS-1	M0		12601		68.005	
			PEP	1,2,3		MS-1	M1		2301		16.505	
			PEP	1,2,3		MS-1	M2		96		5.48	
			PEP	1,2,3		MS-1	M3		1		5.005	

Here, column ``Time`` is left empty intentionally, thus signaling a stationary labeling.

.mflux
------
Starting from this file format, we consider that column names are either similar to already described or are self-explanatory, and we will just give multi-line examples with few possible comments. So, for net flux stationary measurements, we could have: ::

	Id	Comment	Flux	Value	SD
			upt	1.02	0.05


.mmet
-----
For stationary specie concentration, we could have: ::
  
   Id	Comment	Specie	Value			SD
   		Fru6P	0.4263681348074568	0.01
   		GA3P	0.469998855791378	0.01


.cnstr
------
Constraints on fluxes and specie concentrations can look like: ::

	Id	Comment	Kind	Formula			Operator	Value
			NET	Glucupt_1+Glucupt_U	==		1
			NET	edd			>=		0.0001

Column ``Kind`` indicates if a constraint is on net fluxes: ``NET``; on exchange fluxes: ``XCH`` or on specie concentrations: ``MET``. The ``Formula`` content must be a linear function of involved entities. If numeric factors are involved in the formula, they must preceed the variable name, e.g. ``0.632*BM`` and not ``BM*0.632``. Column ``Value`` can have either a float number or a simple Python arithmetic expression which evaluates to a float number, e.g. ``math.sqrt(2)/2`` or ``np.sqrt(2)/2`` (here ``np`` stands for ``numpy``).

.tvar
-----
Types of variables can resemble to ::

	Id	Comment	Name	Kind	Type	Value
			upt	NET	F	1.02369
			emp1	NET	F	0.511098
			emp2	NET	D	
			ppp2	XCH	F	0.778786
			ppp3	XCH	C	0.83932

The ``Type`` can be either

  ``F``
     for "free", requires a float number in ``Value``
  ``D``
     for "dependent" or
  ``C``
     for "constrained", requires a float number in ``Value``.

The ``Kind`` can be either

  ``NET``
     for net fluxes
  ``XCH``
     for exchange fluxes
  ``METAB``
     for specie concentration
     

.opt
----
Options passed to ``influx_si`` can be similar to: ::

	Id	Comment	Name			Value
			dt			1
			nsubdiv_dt		4
			file_labcin		e_coli_msen.txt
			commandArgs		--noscale --TIMEIT --time_order=2 --zc=0 --clowp 1.e-9
			optctrl:nlsic:errx	1.e-3
			optctrl:nlsic:maxit	50
			optctrl:nlsic:btmaxit	16
			optctrl:nlsic:btstart	1
			optctrl:nlsic:btfrac	0.5
			optctrl:nlsic:adaptbt	1
			optctrl:nlsic:monotone	1
			posttreat_R		plot_ilab.R

The meaning of each possible option is described in different sections of this manual.

.vmtf
-----
Variable part of MTF approach can be used to combine constant and variable
parts of experiments to launch a calculation of flux maps in a batch.
E.g. in a set of experiments on the same organism in different biological conditions, ``.miso``, ``.mflux`` can vary from one experiment to another while ``.netw`` and other files can remain the same in the whole experiment set.
In this case, files containing variable sections (here ``.miso`` and ``.mflux``)
can be given in a special file having an extension '.vmtf' while constant
parts will be given in ``--mtf`` or ``--prefix`` options.

``.vmtf`` file is a TSV file with 
columns using the same name as extensions described above: ``netw``, ``linp``, etc.
Each row contains file names that will be used to produce a particular
FTBL file used in calculation.
Thus, each row must have ``ftbl`` column with unique and non empty name. If a file 
type is present both in column names of 'vmtf' and in ``--mtf``/``--prefix`` option 
then the content of 'vmtf' file will take precedence. Empty values 
in ``vmtf`` file are ignored. All file paths in ``vmtf`` file are 
considered relative to the location of ``vmtf`` file itself. If in ``.vmtf``,
a file name is given without extension, it is deduced from column name. Example of ``.vmtf file``: ::

	Id	Comment	miso		mflux		tvar		ftbl
			model_WT_BW_1	model_WT_BW_1	model_WT_BW_1	vmtf_WT_BW_1
			model_WT_BW_2	model_WT_BW_2	model_WT_BW_2	vmtf_WT_BW_2
			model_zwf_1	model_zwf_1	model_zwf_1	vmtf_zwf_1


Example of command line using ``.vmtf``: ::

  --prefix ecoli --mtf variable.vmtf
  
Output format
~~~~~~~~~~~~~

Since v7.0, MTF is also used as output format for most result files. However, compared to input files, some additional columns were added. We will describe them in appropriate subsections.

The whole set of output files go to ``mynetwork_res`` directory if this name is not overwritten with ``--out`` option. Here ``mynetwork`` is just a name example to be adapted to your situation. Let review output files.
By default, the result subdirectory is located in the same directory as input MTF files.

  
.. note::

   All result files are silently overwritten if already exist.
   So take care to copy your results elsewhere if you want to protect them from overwriting.

``mynetwork.log``
-----------------

contains the run-time log messages, in particular,
it contains a report on convergence history during the fitting process.
It can be helpful for identifying potential problems. Some warnings can be written in this file, so user should scrutinize this file for lines starting with ``***Warning:`` to be informed about potential issues.

For example, a warning about structurally non identifiable fluxes can look like:

.. code-block:: text

 ***Warning: inverse of covariance matrix is numerically singular.
 Statistically undefined parameter(s) seems to be:
 f.x.pyk
 For a more complete list, see SD column in '.tvar.sim' result file.

``mynetwork.err``
-----------------

This file contains critical error messages preventing from normal program execution. Normally, it should be empty (0 byte size). If it is not, there is something definitely wrong either in the input files or with the program itself that should be fixed before going on.

Problems can appear in all stages of a software run:

* parsing MTF/FTBL files
* R code writing
* R code execution

	* vector-matrix initialization
	* optimization
	* post-optimization treatment

Most of the error messages are automatically generated by underlying languages Python and R. These messages can appear somewhat cryptic for a user unfamiliar with these languages. But the most critical error messages are edited to be as explicit as possible. For example, a message about badly structurally defined network could be similar to

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

.. note::

  In this error message, we report cumomers whose balance gave a zero row in the cumomer matrix (here ``cit_c:<N>`` cumomers, where <N> is an integer, its binary mask indicates the "1"s in the cumomer definition) as well as a list of fluxes having 0 value. This information could help a user to get insight about a flux whose zero value led to a singular matrix. A workaround for such situation could be setting in the ``.cnstr`` file an inequality constraining a faulty flux to keep a small non zero value. A more radical workaround could be restricting some flux classes (input-output  fluxes with the option ``--cinout=CINOUT`` or even all non-reversible ones with the option ``--clownr=CLOWNR``) to stay out of 0, e.g.:
 
 ``$ influx_s.py --clownr 0.0001 --prefix mynetwork``
 
 Adding such inequalities does not guaranty that cumomer matrix will become invertible, but often it does help.
 It's up to the user to check that an addition of such inequalities does not contradict biological sens of his network.

and so on.

A user should examine carefully any warning/error message and start to fix the problems by the first one in the list (if there are many) and not by the easiest or the most obvious to resolve. After fixing the first problem, rerun ``influx_si`` to see if other problems are still here. Sometimes, an issue can induce several others. So, fixing the first issue could eliminate some others. Repeat this process, till all the troubles are eliminated.

``mynetwork.miso.sim``, ``mynetwork.mflux.sim``, ``mynetwork.mmet.sim``
-----------------------------------------------------------------------

These files contain simulated values and are formatted as input files for isotopic, flux and concentration measurements. Additional columns are:

   ``Residual``
     containing reduced residual values calculated as "(simulated - measured)/SD";

   ``Pvalue``
     result of Z-test on residual value. Lower p-value corresponds to higher (in absolute value) residual.

If multiple starting points were used (e.g. with ``--fseries`` and/or ``--iseries`` options) then result corresponding to each starting point will be numbered with prefix ``.V``, e.g. ``mynetwork.V10.miso.sim`` for tenth staring point.


``mynetwork.tvar.sim``
----------------------

Estimated fluxes and specie concentrations (if any) are all gathered in this file, in column ``Value``. Additional columns are:

  ``SD``
    Estimated Standard Deviation
  ``Struct_identif``
    Contains "yes" if the value on the current row is considered as structurally identifiable ("no" otherwise). We arbitrary set a threshold for "yes" decision at SD <= 10000. Note that a value can be structurally identifiable but still statistically poorly defined.
  ``Low_mc``, ``Up_mc``
    If Monte-Carlo method for sensitivity estimation is used (cf. option ``--sens mc=MC``) then these columns contain lower and upper 95% confidence interval limits. They are estimated at quantile levels 2.5% and 97.5% respectively.
 
In case of multiple starting points, this file is numbered.
 
``mynetwork.stat``
------------------

This file contains results of chi2 statistical test. Its content can look like:

.. code-block:: none

	chi2_value		61.5141593189625
	chi2/df			1.86406543390795
	number_of_measurements	54
	number_of_parameters	21
	degrees_of_freedom	33
	p-value			0.998130284339227
	conclusion		At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD

The field names are self explanatory. In case of multiple starting points, this file is numbered.

Custom files
------------

These files (e.g. ``mynetwork.pdf``) can be produced by user supplied scripts that are executed at the end of ``influx_si`` simulations. For example, we provide a script ``plot_ilab.R`` which can be used to plot label kinetics simulated by ``influx_i`` or ``plot_smeas.R`` for ``influx_s``. One or many of such custom scripts can be given in .opt file, field ``posttreat_R`` (cf. e_coli_i.opt for example)

Old result files
----------------

For users willing to continue usage of old result files, namely ``mynetwork_res.kvh``, we kept them in ``tmp`` subdirectory.

Thus typically, a ``tmp`` subdir contains:

 - ``mynetwork.ftbl``
 - ``mynetwork.pres.csv``
 - ``mynetwork.R``
 - ``mynetwork_res.kvh``
 - ``edge.netflux.mynetwork.attrs``
 - ``edge.xchflux.mynetwork.attrs``

Basic influx_si usage
~~~~~~~~~~~~~~~~~~~~~
``influx_si`` can be run as ::

 $ influx_s.py --prefix mynetwork
 
or ::

 $ influx_i.py --prefix mynetwork

Letters ``s`` and ``i`` stand for "stationary" and "instationary".
We suppose here that a valid MTF file set was created. Moreover, we supposed ``influx_s.py`` and ``influx_i.py`` are in the PATH variable.

In the rest of this manual, we'll use just ``influx_s.py`` as an example if the example is valid for both stationary and instationary contexts. If some usage is valid exclusively for ``influx_i.py``, it will be duly signaled.


In a high throughput context, it can be useful to proceed many MTF set files in parallel. This can be done by giving all variable parts of experiment set in a ``.vmtf`` file, e.g. ::

 $ influx_s.py --prefix mynetwork --mtf variable.vmtf

All files are then proceeded in separate independent processes launched almost simultaneously by a bunch of size equal to the number of available or requested cores (if an option ``--np=NP`` is used). It is an operating system who is in charge to make a distribution of all these processes among all available CPUs and cores.

Sometimes, particular cases need usage of special options of ``influx_si``. The list of available options can be seen by running::

 $ influx_s.py --help

If used with options, ``influx_si`` can be run like ::

 $ influx_s.py [options] --prefix mynetwork

where ``[options]`` is an option list separated by a white character.

.. note::
 Here and throughout this manual, content placed in brackets  ``[...]`` is meant to be an optional part of the command. If the user does wish to type an optional part of a command, the brackets themselves must be omitted.

Each option starts with a double dash ``--`` and can be followed by its argument if applicable. For example, to use BFGS optimization method instead of the default NLSIC algorithm, a user can run::

 $ influx_s.py --meth BFGS --prefix mynetwork

or ::

 $ influx_s.py --meth=BFGS --prefix mynetwork

The option names can be shortened till a non-ambiguous interpretation is possible, e.g., in the previous example, the option could be shortened as ``--me BFGS`` or ``--me=BFGS`` because there is no other option name starting by ``me``. But an option ``--no`` could not be distinguished between ``--noopt`` and ``--noscale``. So at least ``--nos`` (for ``--noscale``) or ``--noo`` (for ``--noopt``) should be provided. There is only one option that does not admit a usage of an equal sign to provide an argument, it is ``--excl_outliers``. Use only a space character to provide an argument to this option when required.

Here after, the available options with their full names are enumerated and detailed.

``influx_si`` command line options
----------------------------------
	--version        show program's version number and exit
	-h, --help       show the help message and exit
	--noopt          no optimization, just use free fluxes as is (after a projection on feasibility domain), to calculate
			dependent fluxes, cumomers, stats and so on
	--noscale        no scaling factors to optimize => all scaling factors are assumed to be 1

			This option can be useful if your measurements are already scaled to sum up to 1 which is often the case of MS data. Then, user saves some free parameters corresponding to scaling factors. This option can become mandatory if user wants to prevent scaling factors to be adjusted by optimization process.
	--meth=METH        method for optimization, one of ``nlsic|BFGS|Nelder-Mead|pso``.
                     Default: ``nlsic``. Multiple occurrences of this
                     option can appear on command line. In this case,
                     specified minimization methods are applied successively,
                     e.g. ``--meth pso --meth nlsic`` means that ``pso`` will be
                     used first, then ``nlsic`` will take over from the point
                     where ``pso`` ends. In case of multiple methods, it is
                     recommended starting with non-gradient methods like ``pso``
                     or ``Nelder-Mead`` and make them follow by gradient based
                     methods like ``nlsic`` or ``BFGS``. If ``pso`` or ``Nelder-Mead``
                     are indeed used as the first method, it is not
                     recommended combining them with ``--zc`` option.
	--fullsys        calculate all cumomer set (not just the reduced one
			necessary to simulate measurements)

			This option influences only post-optimization treatment. The fitting itself is still done with the reduced cumomer set or EMU variables if requested so. See the original paper on ``influx_s`` for more information on the reduced cumomer set.
	--emu            simulate labeling in EMU approach

			This option should not produce a different result in parameter fitting. It is implemented and provided in a hope that on some network the results can be obtained in a shorter time
	--irand          ignore initial approximation for free parameters (free fluxes and specie concentrations) from the FTBL file or from a dedicated file (cf --fseries and --iseries
			option) and use random values drawn uniformly from [0,1]
									 
			It is recommended to use this option in conjunction with "--zc 0" option.
	--sens=SENS      sensitivity method: SENS can be 'mc[=N]', mc stands for
			Monte-Carlo. N is the number of Monte-Carlo simulations.
			Default for N: 10

			The sensitivity information (i.e., the influence of the noise in the data on the estimated parameter variation) based on linearized statistics is always provided. So the user has to use this option only if he wants to compare this linearized information to the Monte-Carlo simulations. Note that the default value 10 for the number of simulations is far from to be sufficient to get reliable statistical estimations. This default option allows only to quickly check that this option is working as expected.
	--cupx=CUPX      upper limit for reverse fluxes. Must be in interval [0, 1]. Default: 0.999
	--cupn=CUPN      upper limit for net fluxes. Default: 1.e3
	--cupp=CUPP      upper limit for specie pool. Default: 1.e5
	--clownr=CLOWNR  lower limit for not reversible free and dependent fluxes.
			Zero value (default) means no lower limit

			A byproduct of this option is that it can drastically reduce  cumomer system sizes. As it ensures that non-reversible fluxes cannot change the sign, revers fluxes can be eliminated from pathways leading to observable cumomers. 
	--cinout=CINOUT  lower limit for input/output free and dependent fluxes.
			Must be non-negative. Default: 0
	--clowp=CLOWP    lower limit for free specie pools. Must be positive. Default 1.e-8
	--np=NP            When integer >= 1, it is a number of parallel threads (on
			Unix) or subprocesses (on Windows) used in Monte-Carlo
			(M-C) simulations or for multiple FTBL inputs. When NP is
			a float number between 0 and 1, it gives a fraction of
			available cores (rounded to closest integer) to be used.
			Without this option or for NP=0, all available cores in a
			given node are used for M-C simulations.
	--ln             Least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient

			Jacobian can become rank deficient if provided data are not sufficient to resolve all free fluxes. It can be useful to determine fluxes that can still be resolved by the available measurements. If the Jacobian does not become rank deficient, this option has no influence on the found solution, neither on the optimization process. But if the Jacobian does become rank deficient, a warning message is printed in the error file even if the optimization process could go to the end.

			.. note:: Use this option with caution, in particular, when used in conjunction with Monte-Carlo simulations. As undetermined fluxes will be given some particular value, this value can be more or less stable from one Monte-Carlo simulation to another. This can create an illusion that a flux is well determined. See the linearized statistics in the result file to decide which fluxes are badly resolved.

			A correct way to deal with badly defined metabolic network is to provide additional data that can help to resolve all the fluxes and/or to optimize input label, not just put ``--ln`` option and cross the fingers.

			.. warning:: In this option, the notion of "least norm" is applied to *increments* during the optimization, not to the final solution. So undetermined fluxes could vary from one run to another if the optimization process is started from different points, while well determined fluxes should keep stable values.
	--sln            Least norm of the solution of linearized problem (and not just of increments) is used when Jacobian is rank deficient
	--tikhreg        Approximate least norm solution is used for increments
			during the non-linear iterations when Jacobian is rank
			deficient
									 
			To obtain an approximate solution, a Tikhonov regularization is used when solving an LSI problem. Only one of the options ``--ln`` and ``--tikhreg`` can be activated in a given run.
	--lim            The same as --ln but with a function limSolve::lsei()
	--zc=ZC          Apply zero crossing strategy with non-negative threshold
			for net fluxes
									 
			This option can accelerate convergence in situations when a net flux has to change its sign during the optimization iterations. Once such flux is identified, it is better to write the corresponding reaction in opposite sens in the FTBL file or to give a starting value with a correct sign to avoid such zero crossing situation.
	--ffguess        Don't use free/dependent flux definitions from FTBL
			file(s). Make an automatic guess.
									 
			The fact that free fluxes are chosen automatically does not allow specifying a starting point for optimization iterations so a random starting point is used (drawn uniformly in [0; 1] interval). An option ``--seed`` can be useful to make the results reproducible.
	--fseries=FSERIES  File name with free parameter values for multiple
			starting points. Default: '' (empty, i.e. only one
			starting point from the FTBL file is used)
										 
			The file must be formatted as plain text file with tab separator. There must be as many columns as starting points and at least as many rows as free parameters assigned in this file. A subset of free parameters can be used in this file. In this case, the rest of parameters take their unique starting values from the FTBL file. The first column must contain the names of free parameters used in this file. If there are extra rows whose names are not in the set of free parameter names, they are simply ignored. The first row must contain the names of starting points. These names can be just numbers from 1 to the number of starting points.
	--iseries=ISERIES  Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default '' (empty, i.e. all provided starting points are used)
										 
			When used with conjunction with ``--fseries``, this option indicates the starting points to use from FSERIES file. But this option can also be used in conjunction with ``--irand`` to generate a required number of random starting points, e.g., ``influx_s.py --irand --iseries 1:10 mynetwork`` will generate and use 10 random starting points.
										 
			For both ``--fseries`` and ``--iseries``, one result file is generated per starting point, e.g., ``mynetwork_res.V1.kvh``, ``mynetwork_res.V2.kvh`` and so on. If starting points comes from a ``--fseries`` then the suffixes ``V1``, ``V2``, ... are replaced by the column names from this file. In addition, a file ``mynetwork.pres.csv`` resuming all estimated parameters and final cost values is written.
	--seed=SEED        Integer (preferably a prime integer) used for
			reproducible random number generating. It makes
			reproducible random starting points (--irand) but also
			Monte-Carlo simulations for sensitivity analysis.
			Default: none, i.e., current system value is used, so
			random drawing will be varying at each run.
	--excl_outliers    This option takes an optional argument, a p-value between
		0 and 1 which is used to filter out measurement outliers.
		The filtering is based on Z statistics calculated on
		reduced residual distribution. Default: 0.01.

		Excluded outliers (if any) and their residual values are reported in the ``mytework.log`` file. Non-available (``NA``) measurements are considered as outliers for any p-value.
		An optional p-value used here does not give a proportion of residuals that will be excluded from the optimization process, but rather a degree of being a valuable measurement. So, closer to zero is the p-value, the fewer data are filtered out. If in contrary, you want to filter out more outliers than with the default p-value, use a value grater than the default value of 0.01, e.g.: ::

				  influx_s.py --excl_outliers 0.02 mynetwork.ftbl

		.. note::

			Don't use an equal sign "=" to give a p-value to this option. Here, only a white space can be used as a separator (as in the example above).
	--nocalc          generate an R code but not execute it.
											
			This option can be useful for parallel execution of the generated R files via ``source()`` function in cluster environment
	--addnoise        Add centered gaussian noise to simulated measurements written to _res.kvh file. SD of this noise is taken from FTBL file
	
			This option can be helpful for generating synthetic FTBL files with realistic simulated measurements (cf. :ref:`How to make FTBL file with synthetic data?<howto>`).
	--mtf MTF             option passed to txt2ftbl. See help there.
	--prefix PREFIX       option passed to txt2ftbl. See help there.
	--eprl EPRL           option passed to txt2ftbl. See help there.
	--force FORCE         option passed to txt2ftbl. See help there.
	-o OUT, --out OUT     output directory. Default: basename of input file without suffix + '_res'. If empty, no output directory is created. In this case log and error messages are directed to standard output/error channels. Non empty OUT can be used when only one input file or MTF set is given.
	--copy_doc         copy documentation directory in the current directory and
                     exit. If ./doc exists, its content is silently overwritten.
	--copy_test        copy test directory in the current directory and exit. If
                     ./test exists, its content is silently owerriten.
	--install_rdep     install R dependencies and exit.

			starting from v5.3, this installation is made in interactive mode. I.e. if the default installation directory (the first one from a list returned by R's ``.libPaths()``) is not writable by the user then ``influx_si`` will try to install the needed packages in the directory defined in R session variable ``R_LIBS_USER``. If this last does not exist, the user is asked for a permission to create it. This behavior is the default one of R's ``install.packages()`` which is used here.
	--TIMEIT          developer option

			Some portions of code are timed, and the results is printed in the log-file. A curious user can use this option without any harm.
	--prof            developer option

			This option provides much more detailed profiling of the execution than ``--TIMEIT`` option. Only developers can be interested in using such information.

All command line options can also be provided in a .opt file. A user can put them in the field ``commandArgs``, e.g.

  .. code-block:: none
  
	Name		Value
	commandArgs	--meth BFGS --sens mc=100 --np 1


If an option is provided both on the command line and in the .opt file, it is the command line that has the priority. In such a way, a user is given an opportunity to overwrite any option at the run time. Nevertheless, there is no way to cancel a flag option (an option without argument) on a command line if it is already set in the .opt file. For example, if ``--fullsys`` flag is set in the .opt file, the full system information will be produced whatever command line options are.

Parallel experiments
~~~~~~~~~~~~~~~~~~~~

.. _prlexp:

Staring from v4.0, ``influx_si`` offers possibility to treat labeling  data from parallel experiments. Parallel experiments for stationary labeling were described in the literature (e.g. cf. "Parallel labeling experiments and metabolic flux analysis: Past, present and future methodologies.", Crown SB, Antoniewicz MR., *Metab Eng.* 2013 Mar;16:21-32. doi: 10.1016/j.ymben.2012.11.010). But for instationary labeling, at the best of our knowledge, ``influx_si`` is the first software offering parallel experiments treatment.

The main interest of parallel experiments is increased precision of flux estimations. This comes at a price of additional work for experiments and data gathering, but the result is often worth the effort. As usual, before doing a real "wet" experiment, it can be useful to run a few  "dry" simulations to see if planned experiments will deliver desired precision.

To deal with parallel experiments, a user have to prepare a series of additional .miso/.linp couples, one per additional experiment. While the "main" .miso/.linp couple can be given in ``--prefix`` or ``--mtf`` options.

Each couple provides input labeling and measured labeling data corresponding to an experiment. 
This file architecture ensures that a network topology, flux, and specie values are common to all experiments, while entry label and measurements on labeled species are proper to each experiment.

When files are ready, you can run ``influx_si`` on them, e.g. in ``test/prl_exp/mtf`` directory run: ::

  $ influx_s.py --pref e_coli_glc1-6n --eprl e_coli_glc2n,e_coli_glc3n,e_coli_glc4n,e_coli_glc5n,e_coli_glc6n
 
In this example, 6 parallel experiments were used, the "main" being described in files ecoli_glc1-6n and 5 additional ones in files going from ``e_coli_glc2n`` to ``e_coli_glc6n``. Note that we used a compact form of ``--eprl`` options as all .miso/.linp couples used canonical suffixes. In that way, giving only a prefix like ``e_coli_glc2n`` was sufficient to find the both corresponding files.

The command can be shortened even more if ``prl_exp`` option is used in ``.opt`` file.
For example if we write in the main set ``e_coli_glc1-6n.opt``: ::

	Name	Value
	prl_exp	e_coli_glc2n,e_coli_glc3n,e_coli_glc4n,e_coli_glc5n,e_coli_glc6n

then the command to run parallel experiments becomes simply: ::

  $ influx_s.py --pref e_coli_glc1-6n

This example set of files correspond to stationary labeling experiments described in "Complete-MFA: Complementary parallel labeling experiments technique for metabolic flux analysis", Robert W. Leighty, Maciek R. Antoniewicz, *Metabolic Engineering* 20 (2013) 49–55 (with only difference that we use simulated and noised data instead of measured ones).

We also provide an example of simulated instationary parallel experiments in the files ``e_coli_GX_prl`` (main files) and ``e_coli_GX_X`` (secondary files) corresponding to simultaneous consumption of glucose and xylose. The network for these simulations was borrowed from "13C metabolic flux analysis of microbial and mammalian systems is enhanced with GC-MS measurements of glycogen and RNA labeling", Christopher P. Long, Jennifer Au, Jacqueline E. Gonzalez, Maciek R. Antoniewicz, Metabolic Engineering 38 (2016) 65–72. The experiment consisted in dynamic labeling by uniformly labeled glucose (main experiment)  and by uniformly labeled xylose (secondary one). Labeling kinetics MS data are given in ``e_coli_GX_MS.miso`` and ``e_coli_GX_X_MS.miso`` files respectively. To play with this example (still in the same directory), you can run: ::
 
 $ influx_i.py e_coli_GX_prl

Note that set of measured specie fragments as well as sampling time points for instationary labeling are not necessary the same for instationary experiments. They do
can differ. That's why a ``.opt`` can be necessary to add to ``.miso/.linp`` couple to form a complete parallel experiment.

It should be made a clear distinction between parallel experiments described in this section and independent experiments calculated in parallel. The main difference is that parallel experiments of this section share all the same flux map (only labeling pattern and data can differ) while independent experiments can have each its own flux map, be they calculated in parallel or sequentially.

Options in .opt file
~~~~~~~~~~~~~~~~~~~~
In this section, we describe different options that can appear in .opt file

.. _optopt:

Optimization options
--------------------
These options can help to tune the convergence process of the NLSIC (or any other chosen algorithm). These options are prefixed with ``optctrl`` which is followed by a particular optimization method name and ended by an option name. For example, ``optctrl:nlsic:errx`` corresponds to the stopping criterion. A corresponding ``.opt`` portion could look like

.. code-block:: none


	Name			Value
	optctrl:nlsic:errx	1.e-3

NLSIC parameters
................

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

PSO parameters
..............

Particle Swarm Optimization (PSO) is a stochastic optimization method. It can help to avoid local minimums but its convergence is very slow. That's why its usage can be particularly useful if combined with a deterministic algorithm like NLSIC. We have implemented PSO method based on the code from CRAN package `pso v1.0.3 <https://cran.r-project.org/package=pso>`_  published in 2012 by Claus Bendtsen (papyrus.bendtsen at gmail.com). The original algorithm was written for box constrained problems. While ``influx_si`` requires a usage of general linear constraints. So we modified the algorithms accordingly. Its parameters with their default values used in ``influx_si`` are following:

        trace=0
                an integer controlling the trace printing. A zero value means no printing
        fnscale=1
                scale factor for minimized function. It is useless in ``influx_si`` context.
        maxit=100
                maximal iteration number to not overcome
        maxf=Inf
                maximal number of a cost function evaluation
        abstol=-Inf
                stopping criterion by absolute tolerance during approximating the searched minimum. This parameter can only be useful if the searched minimal value is known in advance. It is not the case of ``influx_si``
        reltol=0
                stopping criterion by relative change in the found minimal value
        REPORT = 10
                if tracing is enabled, this parameters gives the number of iterations passed between two successive reports
        s=NA,
                swarm size. If NA, it is automatically determined.
        k=3, p=NA, w=1/(2*log(2)), c.p=.5+log(2), c.g=.5+log(2)
                are parameters governing PSO minimization paths. For their significance see the original `pso documentation <https://cran.r-project.org/web/packages/pso/pso.pdf>`_
        d=NA
                domain diameter
        v.max=NA
                maximum allowed velocity
        rand.order=TRUE
                proceed swarm particles in random order or not
        max.restart=Inf
                maximal allowed restarts
        maxit.stagnate=Inf
                maximal successive iterations allowed without a detected decrease in optimization function.
        trace.stats=FALSE
                return or not detailed statistics about the convergence process (not used in ``influx_si``)
        type="SPSO2011",
                which PSO strategy to use. Available options are "SPSO2011" and "SPSO2007". More about this in the original documentation.
        tolineq=1.e-10
                tolerance for violating of linear constraints that can happen mainly due to rounding errors. 

Other optimization methods
..........................

Names and default values for BFGS and Nelder-Mead algorithms can be found in the R help on ``optim()`` function.

.. _growthflux :

Growth flux option
------------------
If present, this option makes ``influx_si`` take into account growth fluxes :math:`-\mu{}M` in the flux balance, where :math:`\mu` is a growth rate and :math:`M` is a concentration of an internal specie M by a unit of biomass. Only species for which this concentration is provided in a .tvar file, contribute to flux balance with a flux :math:`-\mu{}M`.
This flux can be varying or constant during optimization process depending on whether the specie M is part of free parameters to fit or not. Usually, taking into account of this kind of flux does not influence very much on the estimated flux values. So, this option is provided to allow a user to be sure that it is true in his own case.

The option is activated by a field ``include_growth_flux``:

.. code-block:: none

	Name			Value
	include_growth_flux	1

Value 0 cancels the contribution of the growth fluxes to the general flux balance.

Another necessary option is ``mu`` giving the value of `µ`:

.. code-block:: none

	Name	Value
	mu	0.12

Please note that the specie concentrations by a unit of biomass are reported in a file .tvar as:

.. code-block:: none

	Name	Kind	Type	Value
	Fum	METAB	C	2.47158569399681
	Suc	METAB	F	15.8893144279264
	Mal	METAB	F	6.47828321758155
	...	...

Specie names used in this section must be identical to those used in the .netw file and others. "F" is used as an indicator of a varying specie pool. Such varying species are part of fitted parameters. Column "Value" is used as starting value in the optimization process.

One of valuable originality of ``influx_s``, it is a possibility to couple fluxomics and metabolomics in stationary experiments. It can be done because specie pools can influence labeling in two ways:

 * through specie pooling (due to compartmentalization and/or co-elution during chromatography)
 * through growth fluxes.

This last influence is often of low intensity compared to specie transformation fluxes. In literature, it is typically neglected.

Another possibility that was added ``influx_si`` is to provide measured specie concentrations in ``.mmet`` file:

.. code-block:: none

	Specie			Value				SD
	Suc			15.8893144279264*1.e-3/10.7	1.e-2
	Mal			6.47828321758155*1.e-3/10.7	1.e-2
	Rub5P+Rib5P+Xul5P	1.66034545348219*1.e-3/10.7	1.e-2

Like for other measurements, the user has to provide a name, a value, and a standard deviation for each entry. Species listed in this section must be defined in the .netw file and must have type "F" in the ``.tvar``. Numerical values can be simple arithmetic expressions (as in the example above) which are evaluated during file parsing.

When a specie name is given as a sum of species (e.g. ``Rub5P+Rib5P+Xul5P``) it is interpreted as a list of species to be pooled. It is done proportionally to their concentrations. No numerical factor can appear in this sum. At least one of the species from the list must be free (i.e. to have "F" type in the ``.tvar`` file). Otherwise, all species from the list would be considered as having a fixed concentration and providing a measurement for such species would be meaningless.

.. note:: Species having "F" (as "Free") in column "Type" in a .tvar file are treated as fittable parameters. We recall that species in this file are identified as having "METAB" in column "Kind".

An example of an MTF files having specie sections and involving growth fluxes can be found in ``test/mtf/e_coli_growth.*``.

Post treatment option
---------------------

User can specify a name of one or several R scripts that will be automatically executed after non aborted ``influx_si`` run. This option can be useful, for example, for plain saving of calculation environment in a file for later exploring in an interactive R session or for plotting results in a pdf file and so on. A very basic example of such a script is provided in the file ``R/save_all.R`` and its use can be found in the options of ``test/e_coli.opt`` file.

To activate this option, the script names must be provided in the ``.opt`` file, in the field ``posttreat_R`` and separated by ``'; '``, e.g.:

 .. code-block:: text

	Name		Value
	posttreat_R	save_all.R; plot_smeas.R

The script name is interpreted as a relative path to the directory where the original MTF files are located.  If the file is not found there, it is searched for in ``influx_si/R``.
After execution of ``save_all.R``, a file ``e_coli.RData`` is created. This particular example can be used to restore a calculation R environment by launching R and executing::

 > load("e_coli.RData")
 
After that, all variables defined in influx_si at the end of the calculations will be available in the current interactive session.
To be able to launch custom calculations on these variables, the user has to do some preliminary actions. An example of such actions can be found in a file ``preamble.R`` which can be adapted for user's case.

To write his own scripts for post treatments or explore the calculated values in an interactive session, a user have to know some basics about existent variables where all the calculation results and auxiliary information are stored. Here are few of them:

dirw
	is a working directory (where the original FTBL file is)
dirx
	is an executable directory (where influx_s.py is)
baseshort
	is a short name of the input FTBL file (without the suffix ``.ftbl`` neither the directory part of the path)
param
	is the vector of the estimated parameters composed of free fluxes, scaling parameters (if any) and specie concentrations (if any)
labargs
	is an environment with all necessary data for label simulation (e.g. ``v=lab_sim(param, cjac=FALSE, labargs)``) or residual calculation (e.g. ``rres=lab_resid(param, cjac=TRUE, labargs)``)
jx_f
	is a environment regrouping calculated quantities. Here are some of its fields:
	
	``lf`` a list with different fluxes:
	
		``fallnx``
			a vector of all net and exchange fluxes (here, exchange fluxes are mapped on [0; 1[ interval)
		``fwrv``
			a vector of forward and reverse fluxes (reverse fluxes are "as is", i.e. not mapped)
			
	``xsim``
		is an internal state label vector
	``simlab``, ``simfmn`` and ``simpool``
		are vectors of simulated measurements for label, net flux and specie pools respectively (fitting at the best of influx_s' capacity the provided measurements)
	``res``
	 is the reduced residual vector, i.e. (simulated-measured)/SD
	``ures``
	 is the unreduced residual vector, i.e. (simulated-measured)
	``jacobian``
	 as its names indicates, is the Jacobian matrix (d res/d param)
	``udr_dp``
	 is the Jacobian matrix for the unreduced residual vector (d ures/d param)

measurements
 is a list regrouping various measurements and their SD
nb_f
 is a list of various counts, like number of fluxes, parameters to fit, system sizes and so on
nm_list
 is a list of names for various vectors like fluxes, species, label vectors, measurements, inequalities and so on
ui, ci
 are inequality matrix and right-hand side respectively
 
A full list of all available variable and functions can be obtained in an R session by executing::

 > ls()
 
This list of more than 400 items is too long to be fully described here. We hope that the few items succinctly described in this section will be sufficient for basic custom treatments.

An inspiration for your own custom treatments and/or plotting can be found in files ``plot_ilab.R`` and ``plot_smeas.R`` that plot instationary and stationary data respectively in pdf files.

Exclusive ``influx_i`` options
------------------------------
There is only one exclusive option that can be given on a command line:

	--time_order=TIME_ORDER     Time order for ODE solving (1 (default), 2 or 1,2).
		Order 2 is more precise but more time-consuming. The
		value '1,2' makes to start solving the ODE with the first
		order scheme then continues with the order 2.
		
		The scheme order can be important for the precision of flux and concentration estimations. The impact is not direct, but can be very significant. Please note that it can happen that order 1 fits the data with lower cost value function, but it does not mean that the fluxes/concentrations are better estimated.

Other options can occur as fields in a ``.opt`` file.

 ``nsubdiv_dt``
	 integer number of sub-intervals by which every time interval is divided to increase the precision of time resolution.
	 
	 It can happen that the value 1 (default) is sufficient for a satisfactory flux/concentration estimation. User can gradually increase this value (2, 3, ...) in successive ``influx_i`` runs to be sure that better time resolution does not impact parameter estimation. This property is called *grid convergence*. A grid convergence is necessary to overcome the result dependency on the choice of a numerical discretization scheme. A grid convergence can be considered as achieved when changes in estimated parameters provoked by a grid refinement are significantly lower than estimated confidence intervals for these parameters.
 ``dt``
	 a real positive number, defines a time step in a regular grid in absence of values in "Time" column in ``.miso`` file.
	 If a "Time" values are well present for label kinetics, then this parameter has no effect.
	 
	 A regular time grid for label simulations can be useful on preliminary stage when user only elaborates MTF files and wants to see if label simulation are plausible. It can also help to produce simulated measurements (which can be extracted from the ``_res.kvh`` file) for further numerical experiments like studying convergence speed, parameter identifiability, noise impact and so on.
 ``tmax``
	 a real positive number, defines the end of a regular time grid if "Time" is empty or absent in ``.miso``. Parameters ``dt`` and ``tmax`` must be defined in such a way that there will be at least 2 time points greater than 0 in the time grid.
	 
	 If "Time" values are well present in ``.miso`` then this parameter can be used to limit the time grid on which the simulations are done. If the value in ``tmax`` is greater than the maximal time value defined in ``.miso`` file then this parameter has no effect.

	 .. note::
	 
	  It is very important that the values for time, flux, and specie concentrations be expressed in concordant units. It would be meaningless to give time in minutes, fluxes in mM/h/g and concentrations in mM. This will lead to wrong results.
	
	  For example, if the time is expressed in seconds and concentrations in mM/g then fluxes must be expressed in mM/s/g.
	
 ``funlabR``
	 since v5.4, ``influx_i`` is able to simulate label propagation in a metabolically stationary network from a label input varying in time. User can supply R expressions which will calculate fractions of different input label components as functions of time ``t``. Those expressions can be provided in the ``Value`` column of ``.linp`` file but they can need some helper functions. Few of them are defined in a file ``funlab.R`` included in ``influx_si`` but user can need more of them. Thanks to this field, he can define them in a custom R file, who's name can be provided here. There can be given only one file. However, if user-defined functions are spread over several files, they can be included via ``source()`` function called from this one. For this purpose, a predefined variable ``dirw`` pointing to the current working directory can be useful. It is worth mentioning that the file defined in this field will be executed in a particular environment so that variables created during its execution won't affect ``influx_si``'s ones. The path of the file provided in this field is relative to the .netw's one. Example:
     
	 .. code-block:: text

		Name		Value
		funlabR 	e_coli_iv_funlab.R  # the file 'e_coli_iv_funlab.R' is in the same directory that 'e_coli_iv.*' MTF set
            
	 Functions that are available in ``funlab.R`` are following:
     
	  ``ppulses(tp, Tint, Hint=rep_len(c(1., 0.), length(Tint)))``
	    computes a signal in the form of periodic rectangular pulses in time points ``tp``. Each period is composed of one or several intervals whose length in time is given in numeric vector ``Tint`` and heights of signals are given in optional numeric vector ``Hint``. By default, ``Hint`` is a sequence of 1's and 0's. The very first period starts at t=0. Returns a numeric vector of the same length as ``tp``.
	  ``linterp(tp, knots, v)``
	    computes a signal in the form of continuous linear piece-wise functions in time points ``tp``. The limits of linear intervals are defined in numeric vector ``knots`` and values at limits must be given in numeric vector ``v``. All ``tp`` values must lie between ``min(knots)`` and ``max(knots)``. Returns a numeric vector of the same length as ``tp``.
	  ``steplinpath(tp, nu)``
	    computes labeling in a linear pathway of non-reversible reactions under step labeling with fully unlabeled initial state, i.e. starting from 0. The signal is calculated in time points ``tp``. The pathway is defined by numeric vector ``nu`` which represents turn-over rates (i.e. a ratio Flux/Specie_Concentration). All values in ``nu`` must be pairwise different. Returns a numeric matrix of size m x n, where number of rows m=length(nu) and number of columns n=length(tp). So, for example, if a signal only of the third specie is required, the function can be called as ``steplinpath(tp, nu)[3,]``
	  ``steplinpath2(tp, nu, init=double(length(nu)), height=1.)``
	    The same as ``steplinpath()`` above but initial state can be different from 0. It can be defined in numeric vector ``init`` of the same length as ``nu``. The label amplitude can be given in a scalar ``height``.
	  ``ppulseslinpath(tp, nu, Tint, Hint=rep_len(c(1., 0.), length(Tint)), init=double(length(nu)))``
	    computes labeling of linear non-reversible pathway under input composed of periodic rectangular pulses. Labeling for all species is calculated in time points ``tp``. Each period is composed of one or several intervals whose length in time is given in numeric vector ``Tint`` and heights of signals are given in optional numeric vector ``Hint``. By default, ``Hint`` is a sequence of 1's and 0's. Initial label levels can be defined in numeric vector ``init``. Returns a numeric matrix of size m x n, where number of rows m=length(nu) and number of columns n=length(tp).
 ``Value`` column in ``.linp``
	 in this field, user can provide R expression calculating fractions of input label as function of time. While for ``influx_s``, it can only be a constant or python expression evaluating to a scalar at compilation time. Such R expressions can refer to a scalar variable ``t`` representing the current time point and should return a scalar value representing a label level between 0 and 1. The sum of all label levels relative to a given specie must be 1 if a full set of isotopomer is provided by user. If one or many isotopomers are lacking, usual conventions apply for completion to 1 (cf. `MTF format`_).

	 .. code-block:: text
     
		# .linp example
		Specie	Isotopomer	Value
		Gluc_U	111111		ppulses(t, c(T1,T2))  # periodic pulses composed of intervals of length T1 ("labeled") and T2 ("unlabeled")
		Gluc_1	100000		ppulses(t, c(T1,T2))
		
		# .opt example
		Name	Value
		funlabR e_coli_iv_funlab.R  // in this R file variables T1 and T2 are defined
            
	 Input isotopomers that are absent in such ``funlab`` fields are supposed to be 0 all the time (except for above-mentioned conventions).
          
	 If a label level cannot be calculated in one arithmetic operation, several R statements can be placed between curled braces ``{}`` separated by semicolon ``;``. The last operation must be the searched result. In the example above, we could exclude usage of helper file ``e_coli_iv_funlab.R`` by defining ``T1`` and ``T2`` directly in the expressions:

	 .. code-block:: text
     
		# .linp example
		Specie	Isotopomer	Value
		Gluc_U	111111		{T1=2; T2=2; ppulses(t, c(T1,T2))}
		Gluc_1	100000		{T1=2; T2=2; ppulses(t, c(T1,T2))}

Old Result File Fields
~~~~~~~~~~~~~~~~~~~~~~

.. note::

 Starting from v7.0 results are written in MTF format described above. So this section is kept only for legacy reason and is not more necessary for reading.

Generally speaking, the names of the fields in the result KVH file are chosen to be self-explanatory. So there is no so much to say about them. Here, we provide only some key fields and name conventions used in the result file.

At the beginning of the ``mynetwork_res.kvh`` file, some system information is provided. Here, "system" should be taken in two sens: informatics and biological. The information is reported in the fields  ``influx`` and  ``system sizes``. These fields are followed by  ``starting point`` information regrouping ``starting free parameters``,  ``starting cost value``, ``flux system (Afl)`` and ``flux system (bfl)``. Name conventions used in these and other fields are the following:

 net and exchange fluxes
	are prefixed by ``n.`` or ``x.`` respectively
 free, dependent, constrained and variable growth fluxes
	are prefixed by ``f.``, ``d.``, ``c.`` and ``g.`` respectively. So, a complete flux name could look like ``f.n.zwf`` which means `free net ZWF flux`.
	Growth fluxes which depend on constant specie concentrations can be found in constrained fluxes. Constant or variable growth fluxes are postfixed with ``_gr`` (as `growth`) string. For example, a flux ``g.n.Cit_gr`` corresponds to a net growth flux of Citrate specie. The growth fluxes are all set as non-reversible, so all exchange fluxes like ``g.x.M_gr`` or ``c.x.M_gr`` are set to 0.
 scaling factors names
	are formed according to a pattern similar to ``label;Ala;1`` which corresponds to the first group of measurements on Alanine molecule in labeling experiments. Other possible types of experiments are ``peak`` and ``mass``.
 MID vector names
	are looking like ``METAB+N`` where ``METAB`` is specie name and ``N`` goes from 0 to the number of carbon atoms in the considered molecule.
 cumomer names
	follow classical convention ``METAB#pattern_of_x_and_1``, e.g. ``Ala#x1x``
 forward and reverse fluxes
	 are prefixed by ``fwd.`` and ``rev.`` respectively, e.g. ``fwd.zwf`` or ``rev.zwf``
 measurement names
	 have several fields separated by a colon ``:``. For example, ``l:Asp:#xx1x:694`` deciphers like:

		 * ``l`` stands for `labeling` experiment (others possibilities are ``p`` for `peak`, ``m`` for `mass` and ``pm`` for `specie pool`)
		 * ``Asp`` is a specie name
		 * ``#xx1x`` is a measurement identification
		 * ``694`` is a line number in the FTBL file corresponding to this measurement.

The field ``optimization process information`` is the key field presenting the results of an optimization process. The fitted parameters are in the subfield ``par``. Other subfields provide some additional information.

The final cost value is in the field ``final cost``.


The values of vectors derived from free fluxes like dependent fluxes, cumomers, MID and so on are in the corresponding fields whose names can be easily recognized.

Linear stats and Monte-Carlo statistics are presented in their respective fields. The latter field is present only if explicitly requested by user with ``--sens mc=MC`` option. In this kvh section, a term ``rsd`` means "relative standard deviation" (in literature, it is often encountered a synonym CV as Coefficient of Variation). It is calculated as SD/Mean and if expressed in percentage then the formula becomes 100%*SD/Mean.

The field ``jacobian dr_dp (without 1/sd_exp)`` report a Jacobian matrix which is defined as a matrix of partial derivatives :math:`\partial{r}/\partial{p}` where *r* is residual vector (Simulated--Measured) and *p* is a free parameter vector including free fluxes, scaling factors (if any) and free specie pools (if any). Note that in this definition, the residual vector is not yet scaled by standard deviation of measurements. Sometimes, Jacobian is called *sensitivity matrix*, in which case a special care should be brought to the sens of derivation. Often, by sensitivity matrix, we intend a matrix expressing how estimated fluxes are sensitive to variations in the measurement data. Such definition corresponds to generalized inverse of Jacobian and it is reported in the field ``generalized inverse of jacobian dr_dp (without 1/sd_exp)``

Network values for Cytoscape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Several network values formatted for cytoscape are written by ``influx_si`` to their respective files. It can facilitate their visualizing and presentation in graphical mode. All these values can be mapped on various graphical attributes like edge width, node size or color scale of them. These files are written at the end of calculations, so if an error has interrupted this process, no such file will be produced. Take care to don't use an outdated copy of these files.

A file named ``edge.netflux.mynetwork.attrs`` can help to map net flux values on edges of a studied network. A file ``edge.xchflux.mynetwork.attrs`` do the same with exchange fluxes. And finally, ``node.log2pool.mynetwork.attrs`` provides logarithm (base 2) of pool concentrations. They can be mapped on some graphical attribute of network nodes.

See `Additional tools`_ section, `ftbl2xgmml: cytoscape view`_ paragraph to know how to produce files importable in Cytoscape from a given MTF set. User's manual of Cytoscape has necessary information about using visual mapper for teaching how some values like net flux values can be mapped on graphical elements like edge width and so on.

Problematic cases
~~~~~~~~~~~~~~~~~

Obviously, everyone would like be able just run a flux estimation software and simply get results, but unfortunately it does not work in this way every time.
In this section, we review some problematic cases which can be encountered in practice.

Structurally non-identifiable fluxes
------------------------------------

It can happen that collected data are not sufficient to resolve some fluxes in your network. Due to the non-linear nature of the issue, this situation can appear for some set of free flux values and disappear for others, or be persistent for any free flux values. An error is reported to signal such situation, e.g.

.. code-block:: text

 lsi: Rank deficient matrix in least squares
 1 unsolvable variable(s):
 f.n.PPDK        7

and execution is stopped.

Various options are then available for a user facing such situation.

1. Collect more data to resolve lacking fluxes. As a rule of thumb, data must be collected on species which are the nodes of convergence of badly defined fluxes or on species situated downhill of convergence point and preserving labeling pattern. The nature of collected data can also be important. Examples can be constructed where mass data are not sufficient to determine a flux but RMN data can do the job.
 
 Before using real data collection, you can make a "dry run" with ``--noopt`` option and with fictitious or even NA values for intended to collect Isospecie in the .miso file. Thus, we can see if, with these new data, the network becomes well resolved. How? If the error message disappear and SD values in the section ``linear stats`` are not very high then chances are that additionally collected data can help to resolve the fluxes.
 
2. Optimize input label. It can happen that you do collect data on a specie situated in convergence point for undefined fluxes, but incoming fluxes are bringing the same labeling pattern which prevents flux(es) to be resolved. May be changing substrate label can help in this situation. For label optimization you can use a software called IsoDesign, distributed under OpenSource licence and available here http:://metatoul.insa-toulouse.fr/metasys/software/isodes/ (may be you have received ``influx_si`` as part of IsoDesign package, in which case you have it already).
 
 Naturally, this label optimization should be done before doing actual experiments. See IsoDesing tutorial for more details on how to prepare and make such optimization.
 
 If you don't want or don't have a possibility to use a software for label optimization or you think to have an insight on what should be changed in substrate labeling to better define the fluxes, you can still make a try with ``influx_s.py --noopt --prefix mynetwork --mtf new_label.linp`` to see if a new labeling will do the job (here ``new_label.linp`` is an example name for a ``.linp`` file set that you will prepare with new entries. It is important that ``--mtf new_label.linp`` comes after ``--prefix mynetwork`` to take precedence over the old one ``mynetwork.linp``)

3. Use ``--ln`` option. It won't make your fluxes well-defined, it will just continue calculation trying to resolve what can be solved and assigning some particular values (issued from so-called *least norm* solution for rank deficient matrices) to undefined fluxes. You will still have a warning similar to:

 .. code-block:: text

	 lsi_ln: Rank deficient matrix in least squares
	 1 free variable(s):
	 f.n.PPDK        7
	 Least L2-norm solution is provided.
 
 informing you that some flux(es) in the network is(are) still undefined. This option can be helpful if undefined fluxes are without particular interest for the biological question in hand and their actual values can be safely ignored.

4. You can give an arbitrary fixed value to an undefined flux by declaring it as constrained in the ``.tvar`` file (letter ``C`` in the column ``Type`` followed by some value in ``Value`` column).

Badly defined fluxes
--------------------

Also known as *statistically undefined fluxes*, these fluxes have big or even huge SD values. The difference between these fluxes and structurally undefined fluxes is that the badly defined fluxes can become well defined if the noise is reduced or hypothetically eliminated. While the latter will still be undetermined even in the absence of the noise. Despite this difference, all options presented in the previous section are applicable here (all but ``--ln`` which would be without effect here).

An additional measure can be taken which consist in experimental noise reduction. Generally, it can be done by using better protocols, better instruments or simply by increasing the measurement repetition number.

Once again, a use of ``--noopt`` with new hoped SD values in the ``.miso`` file can help to see if these new measurements with better noise characteristics will resolve or not the problem.

Slow convergence
----------------

Slow optimization convergence can manifest by following warnings::

 nlsic: Maximal non linear iteration number is achieved

or/and ::

 nlsic: Maximal backtrack iteration number is achieved
 
Theoretically, user can increase the limit for those two numbers
(``optctrl:nlsic:maxit`` and ``optctrl:nlsic:btmaxit`` respectively in the ``.opt`` file) but generally it is not a good idea. It can help only in very specific situations that we cannot analyze here, as we estimate them low probable.
In all cases, a slow convergence is due to high non-linearity of the solved problem. What can vary from one situation to another, it is the nature of this non-linearity. Depending on this nature, several steps can be undertaken to accelerate optimization:

1. If a non-linearity causing the slow convergence is due to the use of function absolute value :math:`|x|` in the calculation of forward and revers fluxes from net and exchange fluxes, then an option ``--zc=ZC`` (zero crossing) can be very efficient. This non-linearity can become harmful when during optimization a net flux has to change its sign, in other words, it has to cross zero.

 This option splits the convergence process in two parts. First, a minimum is searched for fluxes under additional constraints to keep the same sign during this step. Second, for fluxes that reached zero after the first step, a sign change is imposed, and a second optimization is made with these new constraints.
 If ``--zc`` option is used with an argument 0 (``--zc=0`` or ``--zc 0``), it can happen that fluxes reaching zero produce a singular (non invertible) cumomer balance matrix. In this case, an execution is aborted with an error starting like
 
	.. code-block:: text
	 
		Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3 or constrain some of the fluxes listed below to be non zero
		...
	 
 To avoid such situation, an argument to ``--zc`` must be a small positive number, say ``--zc 0.001``. In this case, positive net fluxes are kept over 0.001 and negative fluxes are kept under -0.001 value. In this manner, an exact zero is avoided.
 
 Another way to avoid problem induced by using module function :math:`|x|` is to add inequality(-ies) imposing sens of reaction in ``.cnstr`` file e.g. ::
	
	Id	Comment	Kind	Formula	Operator	Value
			NET	mae	>=		0

 
 Naturally, in this example, you have to be sure that the reaction catalyzed by malic enzyme (here ``mae``) must go in the sens written in your ``.netw`` file.
 
 You can find potential candidates to impose sens of reaction by examining the flux values in ``mynetwork_res.kvh`` after a slow convergence and looking fluxes whose sign (positive or negative) looks suspicious to you. In our practice, we could observe a dramatic increase in convergence speed and stability just after imposing sens of reaction to a "key" reaction. Obviously, such constraint must be in accordance with biological sens of a studied network and its biological condition.
 
2. A high non-linearity can appear for some particular set of fluxes, especially when they take extreme values. E.g., when exchange fluxes are close to 1 or net fluxes take very high values of order 10² or even 10³ (supposing that the main entry flux is normalized to 1). In such a case, user can low this limits (options ``--cupx=CUPX`` and ``--cupn=CUPN`` respectively) or try to exclude outliers (``--excl_outliers P-VALUE``) as outliers can attract the solution in weird zone of fluxes. In this latter case, the first convergence will continue to be slow and will generate corresponding warnings but the second one (after a possible automatic elimination of outliers) can converge much faster.


Convergence aborted
-------------------
This situation is signaled by an error message::

 nlsic: LSI returned not descending direction

This problem can occur for badly defined network, which are very sensitive to truncation errors. The effect of such errors can become comparable to the effect of the increment step during optimization. It means that we cannot decrease the norm of residual vector under the values resulting from rounding errors.
If it happens for relatively small increments, then the results of convergence are still exploitable. If not, there is no so many actions that user could undertake except to make his system better defined as described in previous sections.

.. note:: By default, we use a very small value for increment norm as stopping criterion (:math:`10^{-5}`). It can be considered as very drastic criterion and can be relaxed to :math:`10^{-3}` or :math:`10^{-2}` depending on required precision for a problem in hand (to do that, use an option ``optctrl:nlsic:errx`` in the ``.opt`` file). 

Additional tools
~~~~~~~~~~~~~~~~

Tools described in this section are not strictly necessary for running ``influx_si`` and calculating the fluxes. But in some cases, they can facilitate the task of tracking and solving potential problems in FTBL preparation and usage.

Most of the utilities produce an output written on standard output or in a file whose name is derived from the input file name. This latter situation is signaled with a phrase "The output redirection is optional" and in the usage examples the output redirection is taken in square brackets ``[> output.txt]`` which obviously should be omitted if an actual redirection is required. Such behavior is particularly useful for drag-and-drop usage.

ftbl2mtf: conversion of FTBL to MTF format
------------------------------------------
For old ``influx_si`` users having their projects in FTBL format, this utility can be an invaluable helper for making the transition to the new MTF format. Here is the help message, which can be seen with ``ftbl2mtf -h``

      .. code-block:: text
      
	usage: ftbl2mtf [-h] [-i] [-f] [-o OUT] ftbl

	Parse ftbl file from first parameter or from stdin (if input file is '-')
	and write a series of mtf (multiple TSV files).
	The file stem ('network' in 'network.ftbl') is used as file name basis
	for produced files, e.g. 'network.miso'. Parameter --out can be used to change it.
	If out path includes non existing directories, they are automatically created.
	Caution! If an existing output file starts with a comment
	"# Created by 'ftbl2mft ..."
	or is empty, it is silently overwritten.
	Otherwise, the writing is aborted with a warning. Other files may continue to be created.
	To force the overwriting, use '--force'.

	Output files will have following extensions/meanings:
	
	 .netw: stoichiometric equations and label transitions in the biochemical network;
	 .linp: label input;
	 .miso: isotopic measurements (MS, label, peak);
	 .mflux: flux measurements;
	 .mmet: biochemical specie concentration measurements;
	 .tvar: flux/specie types partition (free, dependent, constrained) and starting values;
	 .cnstr: constraints (equalities, inequalities for both fluxes and concentrations);
	 .opt: options.

	Copyright 2022 INRAE, INSA, CNRS
	Author: Serguei Sokol (sokol [at] insa-toulouse [dot] fr)

	positional arguments:
	  ftbl               input file to be converted to MTF

	optional arguments:
	  -h, --help         show this help message and exit
	  -i, --inst         activate instationary mode
	  -f, --force        force overwriting of result files
	  -o OUT, --out OUT  path prefix for result files



txt2ftbl: conversion of MTF format to FTBL format
-------------------------------------------------
This tool is implicitly used by ``influx_si`` to convert MTF to FTBL format. Users desiring to play with format conversion or to produce FTBL file to be used with `Additional tools`_ can use it explicitly. Here is its help message:

    .. code-block:: text
	
	usage: txt2ftbl [-h] [--mtf MTF] [--prefix PREFIX] [--eprl EPRL] [--inst] [--force] [netw]

	transform a series of TXT and TSV files into FTBL file.

	Copyright 2021, INRAE, INSA, CNRS
	Author: Serguei Sokol (sokol at insa-toulouse dot fr)
	License: Gnu Public License (GPL) v2 http://www.gnu.org/licenses/gpl.html

	positional arguments:
	  netw             
			   If 'netw' file is not given in any option (neither --mtf nor --prefix), it 
			   can be given as the only argument NETW, e.g.
			     txt2ftbl ecoli.txt
			   or
			     txt2ftbl --mtf ms_nmr_data.miso,glucose.linp ecoli.txt
			   If 'netw' file name is given both in any option and as an argument, it 
			   is the argument value that will take precedence.

	optional arguments:
	  -h, --help       show this help message and exit
	  --mtf MTF        MTF is a coma separated list of files with following extensions/meanings:
			    netw: a text file with stoichiometric reactions and label transitions (one per line)
			       Comments starts with '#' but those starting with '###' introduce 
			       pathways which are numbered as well as reactions in them. Reaction 
			       name can precede the reaction itself and is separated by ":" If no 
			       explicit name is given, reactions in FTBL file will be named 
			       according a pattern 'rX.Y' where X is pathway number and Y is 
			       reaction number in the pathway. But it is highly recommended to 
			       give explicit names to reactions.
			       Symbols "+", "(", ")" and ":" are not allowed in metabolite neither reaction names
			       Example of reaction and label transition:
				  edd: Gnt6P (ABCDEF) -> Pyr (ABC) + GA3P (DEF)
			       Non reversible reactions are signaled with '->' (as in the example above).
			       A sign '<->' can be used for reversible reactions.
			       If 'netw' name is equal to '-', then its content is read from standard input, e.g.
				 '--mtf netw=-'
			    linp: label inputs (starting from this extensions, TSV files are assumed)
			    miso: isotopic measurements (NMR (label, peak) and MS)
			    mflux: flux measurements
			    mmet: metabolite concentration measurements
			    tvar: type of variables (NET or XCH , free or dependent, starting values, ...)
			    cnstr: equality and inequality constraints on fluxes and concentrations
			    opt: options
			    ftbl: name of output FTBL file. If not given, it will be equal to 'netw'
			      stem with '.ftbl' extension. If it is equal to '-', then the result 
			      will be written to standard output, e.g.
				'txt2ftbl --mtf ftbl=-,ecoli.netw'
			      Intermediate directories in ftbl path are silently created if non existent.
			    vmtf: variable part of mtf approach.
			      If a series of FTBL files has to be generated partially with 
			      information common to all files (constant part) and partially with 
			      sections proper to each FTBL (variable part) then files containing 
			      variable sections (e.g. 'miso') can be given in a special file 
			      having an extension (or prefix, cf. hereafter) 'vmtf'. In such a 
			      way, each FTBL file will be produced from combination of MTF files 
			      given directly in this option (constant part) and files given on a 
			      corresponding row of 'vmtf' file. vmtf file is a TSV file with 
			      columns using the same names: 'netw', 'linp', etc. Each row contains 
			      file names that will be used to produce an FTBL file. Thus each row 
			      must have 'ftbl' column with unique and non empty name. When 'vmtf' 
			      is used, 'ftbl' cannot be present on the command line. If a file 
			      type is present both in column names of 'vmtf' and in '--mtf' option 
			      then the content of 'vmtf' file will take precedence. Empty values 
			      in 'vmtf' file are ignored. All file paths in 'vmtf' file are 
			      considered relative to the location of 'vmtf' file itself.
			   Only first 3 files are necessary to obtain a workable FTBL file, others 
			   are optional.
			   Example: 'txt2ftbl --mtf ecoli.netw,glu08C1_02U.linp,cond1.miso,cond1.mflux'
			   NB: no space is allowed around comas. If a file path has a spaces in 
			   its name, it must be enclosed into quotes or double quotes.
			   If an entry file cannot be renamed to have some of these extensions, then
			   they can be used as prefixes followed by a '=' sign, e.g.
			     'txt2ftbl --mtf netw=ecoli.txt,linp=glu08C1_02U.tsv,cond1.miso,cond1.mflux'
			   As you can see from this example, both naming schemes can be mixed.
			   If for some reason, the same type of file is indicated several times
			   (no matter with extension or prefix), the last occurrence supersedes
			   all precedent ones.
	  --prefix PREFIX  If all input files have the same name pattern and are different only 
			   in extensions then the pattern can be given as PREFIX, e.g.
			     '--prefix somedir/ecoli'
			   Then in 'somedir', we suppose to have 'ecoli.netw', 'ecoli.linp' and 
			   other input files having names starting with 'ecoli' and ending with 
			   corresponding extensions.
			   NB. If some file is given in more than one option: '--prefix' and/or 
			   '--mtf' then the last occurrence overrides precedent ones.
	  --eprl EPRL      Parallel experiments can be given with this option. It must
			   introduce a couple of linp/miso files and optional auxiliary ftbl name. These files
			   correspond to a given parallel experiment. This option can be repeated as
			   many times as there are additional parallel experiments, e.g.
			     'txt2ftbl --mtf ec.netw,glc6.linp,glc6.miso --eprl glc1.linp,glc1.miso --eprl glc4.linp,glc4.miso'
			   This command will produce a main FTBL file 'ec.ftbl' including all necessary
			   sections (NETWORK, etc.) but also two auxiliary FTBL files: 'glc1.ftbl' and
			   'glc4.ftbl' having only label input/measurement sections. They will correspond
			   to 2 additional parallel experiments. If ftbl file is not given in --eprl
			   option, the name of miso file will be used for it. If intermediate
			   directories in ftbl path are non existent they will be silently created.
			   Auxiliary ftbl names will be put in 'OPTIONS/prl_exp' field on the main ftbl file.
			   These names will be written there in a form relative to the main ftbl.
			   To shorten the writings, it is possible to indicate only one of two .miso/.linp files.
			   The other one will be guessed if it has canonical extension. If extension is omitted then .miso and .linp files are searched with these extensions. In this case, several parallel experiments can be given with one --eprl option. So that above example can be shorten to:
			     'txt2ftbl --mtf ec.netw,glc6.linp,glc6.miso --eprl glc1,glc4'
	  --inst           Prepare FTBL for instationary case. File 'netw' is supposed to have 
			   column 'Time' non empty. Isotopic kinetic data will be written to a TSV 
			   file with 'ikin' extension. Its name will be the same as in FTBL file, 
			   and FTBL field 'OPTIONS/file_labcin' will contain 'ikin' file name.
	  --force          Overwrite an existent result file not produced by this script.
			   NB. If a result file exists and is actually produced by this script, 
			   then it is silently overwritten even without this option. The script 
			   detects if it was the creator of a file by searching for a string "// 
			   Created by 'txt2ftbl" at the first line of the file. By removing or 
			   editing this comment, user can protect a file from a silent 
			   overwriting.

ftbl2xgmml: Cytoscape view
--------------------------

Once a valid MTF set (or FTBL file) is generated, a user can visualize a graph representing his metabolic network in Cytoscape_ program. To produce necessary graph files, user can run::

 $ ftbl2xgmml.py --prefix mynetwork [> mynetwotk.xgmml]
 
or::

 $ ftbl2xgmml.py mynetwork[.ftbl] [> mynetwotk.xgmml]

or drag and drop ``mynetwork.ftbl`` icon on ``ftbl2xgmml.py`` icon.

The output redirection is optional.

This will produce a file in the XGMML format ``mynetwork.xgmml`` in the directory of ``mynetwork.netw`` or ``mynetwork.ftbl``:

Once a generated file ``mynetwork.xgmml`` is imported in cytoscape, a user can use one of automatic cytoscape layouts (e.g. 'Prefuse Force Directed Layout') or edit node's disposition in the graph by hand.
For those who use `CySBML <http://apps.cytoscape.org/apps/cysbml>`_ plugin, a saving of a particular layout in a file can be practical for later applying it to a new network.

Graphical conventions used in the generated XGMML are the following:

* specie are presented as rounded square nodes;
* simple (one to one) reaction are represented by simple edges;
* condensing and/or splitting reactions are represented by edges converging and/or diverging from an additional almost invisible node having a label with the reaction name;
* all nodes and edges have tool tips, i.e., when a pointer is put over, their name (specie or reaction) appears in a tiny pop-up window;
* non-reversible reactions are represented by a single solid line, have an arrow on the target end (i.e., produced specie) and nothing on the source end (i.e., consumed specie);
* reversible reactions are represented by a double parallel line and have a solid circle on the source end;
* color code for arrows:

	* green for free net flux;
	* blue for dependent net flux;
	* black for constrained net flux;

* color code for solid circles:

	* green for free exchange flux;
	* blue for dependent exchange flux;
	* black for constrained exchange flux.

ftbl2netan: MTF/FTBL parsing
----------------------------

To see how MTF/FTBL files are parsed and what the parsing module "understands" in the network, a following command can be run::

 $ ftbl2netan.py --prefix mynetwork [> mynetwork.netan]
 
or::

 $ ftbl2netan.py mynetwork[.ftbl] [> mynetwork.netan]

The output redirection is optional.

A user can examine ``mynetwork.netan`` in vkvh_ software or in a plain text editor (not like Word) or in spreadsheet software. It has an hierarchical structure, the fields are separated by tabulations and the field values are Python objects converted to strings.

ftbl2cumoAb: human-readable equations
-------------------------------------

Sometimes, it can be helpful to examine visually the equations used by ``influx_si``. These equations can be produced in human-readable form by running::

 $ ftbl2cumoAb.py -r --prefix mynetwork [> mynetwork.sys]
 
or::

 $ ftbl2cumoAb.py -r mynetwork[.ftbl] [> mynetwork.sys]

or::

 $ ftbl2cumoAb.py --emu mynetwork[.ftbl] [> mynetwork.sys]
 
The output redirection is optional.

The result file ``mynetwork.sys`` will contain systems of stoichiometric and cumomer balance equations as well as a symbolic inversion of stoichiometric matrix. I.e., dependent fluxes are represented as a linear combination of free and constrained fluxes and an optional constant value. In the examples above, the option ``-r`` stands for "reduced cumomer set" and ``--emu`` stands for "generate EMU framework equations". In this latter case, only isotopologues of mass+0 in each EMU are reported in ``mynetwork.sys`` file. For other mass weights, the equations does not change and the right-hand side term could get longer for condensation reactions but involves the same EMUs as in mass+0 weight.

If a full cumomer set has to be examined, just omit all options. Keep in mind that on real-world networks this can produce more than a thousand equations by cumomer weight, which could hardly be qualified as *human*-readable form. So use it with caution.

For the sake of brevity, cumomer names are encoded in decimal integer form. For example, a cumomer ``Metab#xx1x`` will be referred as ``Metab:2`` because a binary number ``0010`` corresponds to a decimal number ``2``. The binary mask ``0010`` is obtained from the cumomer mask ``xx1x`` by a plain replacement of every ``x`` by ``0``.

For a given cumomer weight, the equations are sorted alphabetically.

expa2ftbl: non-carbon carrying fluxes
-------------------------------------

Deprecated since v6.0. Such kind of fluxes can be directly incorporated in .netw file.

Some reactions of carbon metabolism require cofactor usage like ATP/ADP and some others. A mass balance on cofactors can produce additional useful constraints on the stoichiometric system. Since the version 2.8, such mass balance equation on non carbon carrying species can be put in ``EQUATION`` section of FTBL file. A utility ``expa2ftbl.R`` can be helpful for this purpose if a user has already a full set of reactions in `expa <http://gcrg.ucsd.edu/Downloads/ExtremePathwayAnalysis>`_ format.
To extract additional equation from an expa file, ``expa2ftbl.R`` can be used as::

 $ R --vanilla --slave --args file.expa < expa2ftbl.R > file.ftbl_eq

Then an information for the generated ``file.ftbl_eq`` has to be manually copy/pasted to a corresponding FTBL file.

Note that ``expa2ftbl.R`` uses a Unix command ``grep`` and another utility described here above ``ftbl2netan.py``.

res2ftbl_meas: simulated data
-----------------------------

.. _ex_sim1:

.. note:: Deprecated since v7.0. User can use ``.miso.sim``, ``.mflux.sim`` and ``.mmet.sim`` as measurement files in a new calculation. E.g.

  .. code-block:: none
  
    influx_s --pref mynetwork --mtf \
    miso=mynetwork_res/mynetwork.miso.sim,\
    mflux=mynetwork_res/mynetwork.mflux.sim \
    -o tmp

  here, new output directory ``tmp`` is used to avoid an overwriting of previous results.

During preparation of a study, one of the questions that biologist can ask is "Will the intended collected data be sufficient for flux resolution in a given network?"
Some clue can be obtained by making "dry runs" of ``influx_si`` with ``--noopt`` (i.e. no optimization) option. User can prepare an FTBL file with a given network and supposed data to be collected. At first, the measurement values can be replaced by NAs while the SD values for measurements must be given in realistic manner. After running::

 $ influx_s.py --noopt mynetwork

a utility ``res2ftbl_meas.py`` can be practical for preparing FTBL files with obtained simulated measurements::

 $ res2ftbl_meas.py res2ftbl_meas.py mynetwork_res[.kvh] > mynetwork.ftbl_meas

(here ``.kvh`` suffix is optional). The information from the generated file ``mynetwork.ftbl_meas`` has to be manually copy/pasted into the corresponding FTBL file.
Getting an ftbl file with real values instead of NAs in measurement sections gives an opportunity to explore optimization behavior near a simulated point like convergence speed and/or convergence stability to cite few of them.

ffres2ftbl: import free fluxes
------------------------------

.. _ex_sim2:

.. note:: Deprecated since v7.0. User can use ``.tvar.sim`` as starting values in a new calculation. E.g.

  .. code-block:: none
  
    influx_s --pref mynetwork --mtf \
    tvar=mynetwork_res/mynetwork.tvar.sim \
    -o tmp

  here, new output directory ``tmp`` is used to avoid an overwriting of previous results.

This utility imports free flux values and specie concentrations (if any) from a result file _res.kvh and inject them into an FTBL file. Usage::

 $ ffres2ftbl.sh mynetwork_res.kvh [base.ftbl] > new.ftbl

If an optional argument ``base.ftbl`` is omitted, then the free flux values are injected into an FTBL file corresponding to the _res.kvh file (here ``mynetwork.ftbl``). This script can be used on a Unix (e.g., Linux, MacOS) or on a cygwin (Unix tools on Windows) platform. It makes use of another utility written in python ``ff2ftbl.py``

ftbl2kvh: check MTF/FTBL parsing
--------------------------------

This utility simply parses MTF/FTBL files and write what was "understood" to a kvh file. No network analysis occurs here unlike in ``ftbl2netan`` utility. Usage::

 $ ftbl2kvh.py --prefix mynetwork [> mynetwork.kvh]
 
or::

 $ ftbl2kvh.py mynetwork[.ftbl] [> mynetwork.kvh]

The output redirection is optional.
The resulting ``mynetwork.kvh`` is tab-separated plain text file. So it can be explored in vkvh_ software, plain text editor or spreadsheet software.

ftbl2metxml: prepare MetExplore_ visualization
----------------------------------------------

Convert a MTF sets (or FTBL files) to an xml files suitable for visualization on MetExplore_ site. If a result kvh file ``mynetwork_res.kvh`` is present, it will be parsed to extract flux values corresponding to the last ``influx_si`` run and put them in ``mynetwok_net.txt``, ``mynetwork_fwd.txt`` and ``mynetwork_rev.txt``. As their names indicate, they will contain net, forward and revers flux values respectively.

IsoDesign: optimizing input label
---------------------------------

One of the means to increase a flux resolution can be an optimization of input label composition. A utility ``IsoDesing`` solving this problem was developed by Pierre Millard. It is not part of ``influx_si`` distribution and can be downloaded at http://metatoul.insa-toulouse.fr/metasys/software/isodes/. In a nutshell, it works by scanning all possible input label compositions with a defined step, running ``influx_si`` on each of them. Then, it collects the SD information on all fluxes for all label compositions and finally selects an input label composition optimal in some sens (according to a criterion chosen by a user).

