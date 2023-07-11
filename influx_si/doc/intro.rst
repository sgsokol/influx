
.. _Python: https://www.python.org/

.. _R: https://www.r-project.org/

============
Introduction
============

``influx_s`` and ``influx_i`` are programs written in Python_, R_ and C++ and designed for estimation of flux and chemical specie concentrations based on labeling data using stable isotopes (essentially ¹³C but combination of multiple isotopes like ²H, ¹³C, ¹⁵N, ... is also possible). ``influx_s`` works with stationary data while the ``influx_i`` is able to simulate instationary labeling (hence the ``_s`` and ``_i`` in the names). Both work in metabolically stationary context. The whole project is referred as ``influx_si``. Note also that the term ``influx_si`` is used in contexts where a proposition holds for both ``influx_s`` and ``influx_i``.

``influx_si``
-------------

Flux and metabolite concentration values are obtained
as a result of a fitting between simulated labeling data and the data measured by MS or NMR techniques. In this documentation, the terms `fitting` and `optimization`
are used as synonyms.

``influx_s``
------------

For the theory behind flux calculations in stationary labeling context, see the following papers:

Wiechert, W., Möllney, M., Isermann, N., Wurzel, M., and de Graaf, A. A. (1999).
Bidirectional reaction steps in metabolic networks: III. Explicit solution and analysis
of isotopomer labeling systems. *Biotechnol Bioeng,* 66(2), 69-85.

Antoniewicz, M. R., Kelleher, J. K., and Stephanopoulos, G. (2007). Elementary
metabolite units (EMU): a novel framework for modeling isotopic distributions.
*Metab Eng*, 9(1), 68-86.

Sokol, S., Millard, P., and Portais, J-C. (2012). 
influx_s: increasing numerical stability and precision for
metabolic flux analysis in isotope labeling experiment.
*Bioinformatics*, 2012, 28, 687-693

The main additional value to flux calculation of ``influx_si`` compared to other publicly
available software (`13CFlux <https://www.13cflux.net>`_,
`OpenFlux <http://openflux.sourceforge.net/>`_, `INCA <http://
mfa.vueinnovations.com>`_, ...) is the usage of NLSIC algorithm
for fitting purposes. This algorithm provides:

 - more reliable convergence which results in better numerical precision, i.e. even started from random initial points, it converges to the same solution if no local minima are present. So the spread of final solutions is close to zero.
 - better accuracy, i.e., the found numerical solution lies closer to the theoretical solution than solutions provided by concurrent minimization algorithms. Thus, ``influx_s`` provides better numerical accuracy.

For more details, see the paper on ``influx_s`` cited above.

Moreover, ``influx_s`` provides:

 - both cumomer and EMU frameworks for describing label distribution in the metabolites;
 - parallel experiment (i.e. same flux/concentration map but different labeling strategies) both in stationary and instationary modes;
 - estimation of specie concentration in particular in stationary contexts (since v2.0. A methodology behind metabolite concentration evaluation is not yet published at the moment of this writing.); 
 - a possibility to deal with metabolite pool confusion appearing either in compartmentation or in coelution;
 - taking into account non-carbon carrying fluxes like the balances of ADP/ATP, H2O, energy, electrons and so on;
 - an optional automatic choice of free fluxes;
 - optional equality and inequality constraint on fluxes and metabolite concentrations;
 - short time execution and design for many core computers. So it facilitates high throughput flux calculations in parallel way;
 - a 'least norm' option that, in presence of structurally non identifiable fluxes, still allows to estimate some of fluxes (those remained identifiable);
 - a chi2 statistical test 'goodness of fit'
 - an optional automatic elimination of outliers;
 - a command line interface letting an easy integration in automatic processing chains as well as many others features and options;
 - a possible scripting of post-treatment or graphic generating tasks;
 - multi-platform support. It runs everywhere R_ and Python_ run, i.e. on Linux, Windows, MacOS and other Unix variants.

``influx_i``
------------

Instationary labeling is the domain of ``influx_i``.
The theory of instationary labeling was developed, for example in

Katharina Nöh, Wolfgang Wiechert (2006)
Experimental Design Principles for Isotopically Instationary 13C Labeling Experiments
*Biotechnology and Bioengineering*, 94(2), 234-251

Sokol S, Portais J-C (2015)
Theoretical Basis for Dynamic Label Propagation in Stationary Metabolic Networks under Step and Periodic Inputs.
*PLoS ONE* 10(12): e0144652. doi:10.1371/journal.pone.0144652

As ``influx_i`` capitalizes on ``influx_s`` development and shares a big part of code, ``influx_i`` presents the same advantages as listed in the previous section. It uses the same input/output file formats for network, measurements definitions and simulated results. It includes all options available for ``influx_s``. Instationary labeling data can be supplied by giving non empty values in ``Time`` column in input file ``.miso`` thus making a shift from stationary to instationary calculations as simple as possible.
Some of advantages of ``influx_i`` over other software coping with instationary labeling data are:

 - fast calculations (e.g. on our Intel Xeon 2.50GHz workstation, ``e_coli_i`` case runs in 17s while the most important part devoted to optimization takes as low as 10s);
 - parallel experiment treatment;
 - available choice between first and second order time schemes for ODE (ordinary differential equations) resolution;
 - unconditional stability during ODE solving.
 
Documentation organization
--------------------------

Changes brought to every new version and bug fixes are resumed at the beginning of the next chapter :doc:`Change Log<changelog>`.

The rest of the documentation is organized as follows. :doc:`Installation <install>` chapter provides brief instructions for software installation. :doc:`Quick start <quick>` chapter gives an opportunity to a user to quickly start and evaluate the software and to see if it corresponds to what he is looking for. A more detailed but still short :doc:`User's manual <manual>` precedes a :doc:`Programmer's documentation <progdoc>`. The latter chapter can be safely skipped by a user not interested in developing new features or fixing some problems in ``influx_si``. A small collection of :doc:`How to... <howto>` and :doc:`Troubleshooting <trouble>` notice conclude the documentation.

Licensing
---------

The original version of ``influx_si`` software was developed in the MetaSys team in the LISBP (TBI since 2018), Toulouse, FRANCE.

The software is licensed under the GNU Public License, Version
2.0 or higher at your convenience (the "License"); you may not use this software and documentation except in compliance with the License.

A file ``influx_si/R/psoptim_ic.R`` is based on the code from CRAN package `pso v1.0.3 <https://cran.r-project.org/package=pso>`_  published in 2012 by Claus Bendtsen (papyrus.bendtsen at gmail.com). The original code is licensed under LGPL-3 terms so our modifications are licensed under the `same terms <https://www.gnu.org/licenses/lgpl-3.0.en.html>`_ .

If you publish results obtained with ``influx_s`` you have to cite the original paper in Bioinformatics 2012 (cf. above). A paper describing ``influx_i`` is yet to publish.

You may obtain a copy of the License :doc:`here <license>` or at

https://www.gnu.org/licenses/old-licenses/gpl-2.0.html

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.


Software and documentation author:

  Serguei SOKOL, INRAE, France <sokol [at] insa-toulouse.fr>

Copyright 2011-2023, INRAE/CNRS/INSA
