============
Introduction
============

``influx_s`` is a software designed for flux calculation based on
labeling data using 13C isotop. The fluxes are calculated in
metabolically and isotopically stationary mode. Their values are obtained
as a result of a fitting between simulated labeling data and the data measured
by MS or NMR techniques. In this documentation the terms `fitting` and `optimization`
are used as synonyms. For the theory behind these calculations see the following papers:

Wiechert, W., Möllney, M., Isermann, N., Wurzel, M., and de Graaf, A. A. (1999).
Bidirectional reaction steps in metabolic networks: III. Explicit solution and analysis
of isotopomer labeling systems. Biotechnol Bioeng, 66(2), 69-85.

Antoniewicz, M. R., Kelleher, J. K., and Stephanopoulos, G. (2007). Elementary
metabolite units (EMU): a novel framework for modeling isotopic distributions.
Metab Eng, 9(1), 68-86.

Sokol, S., Millard, P., and Portais, J-C. (2011). 
influx_s: increasing numerical stability and precision for
metabolic flux analysis in isotope labeling experiment. (Submitted)

The main additional value of this program compared to other publicly
available software (`13CFlux <https://www.13cflux.net>`_, `OpenFlux <http://openflux.sourceforge.net/>`_, ...) is the usage of NLSIC algorithm
for fitting purposes. This algorithm provides:

 - more reliable convergence (even started from random initial points it converges to the same solution if no local minima are present)
 - better accuracy (the found solution has better numerical precision due to a particular stopping criterion).
For more details, see the paper on ``influx_s`` cited here before.

This documentation is organized as follows. :doc:`Installation <install>`
chapter provides brief instructions for software installation. :doc:`Quick start
<quick>` chapter gives an opportunity to a user in a hurry to quickly start and
evaluate the software and to see if it corresponds to what he is looking for. A
more detailed but short :doc:`User's manual <manual>` precedes a
:doc:`Programmer's documentation <progdoc>`. The latter chapter can be safely
skipped by a user not having in a project to develop new features or to fix some
problems in ``influx_s``. A small collection of :doc:`How to... <howto>` and
:doc:`Troubleshooting <trouble>` notice are concluding the documentation.

Licensing
---------

The original version of ``influx_s`` software was developed in the MetaSys team
in the LISBP, Toulouse, FRANCE.

The software is licensed under the Educational Community License, Version
2.0 (the "License"); you may not use this software and documentation except in compliance with the
License. You may obtain a copy of the License :doc:`here <license>` or at

http://www.opensource.org/licenses/ECL-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


Software and documentation author:

  Serguei SOKOL, INRA, France <sokol [at] insa-toulouse.fr>

Copyright 2011, INRA, France
