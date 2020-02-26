
.. _trouble:

===============
Troubleshooting
===============

The software is provided "AS IS" without warranty of any kind explicit or implicit.

If you have some issue with ``influx_si`` you can try the following steps for solving them:

 * partners and clients of MetaToul-RéseauxMétaboliques can benefit from advices of their dedicated stuff. If you are not a partner/client but would like to become one, for example to get help for your label experiment design and/or realization, you can contact platform MetaToul (cf. :doc:`Consulting and more <consulting>`);
 
 * you can search for similar problem discussion in the forum https://groups.google.com/forum/#!forum/influx_si. If you don't find your answer and wish to ask a new question you'll have to subscribe to this group.
 
 * if you think that you face a bug, try the latest version of the software to see if this bug was already fixed. If it is still present, you can report it on https://github.com/sgsokol/influx/issues. Please note that we can't guarantee that any particular bug can be fixed in any particular release or can be fixed at all. It is possible, that we ask you to send us (in a private email not in influx_si@googlegroups.com) an ftbl file on which an error occur. It will be used only for purposes of bug reproducing and its identification. The received ftbl file will not be transmitted to any third party.

 * if you have a problem with FTBL editing, you can read the documentation from `13CFlux <https://www.13cflux.net>`_ and/or interpret error messages generated during FTBL parsing.

 * if you have some difficulties in choosing free fluxes, define all not constrained fluxes as dependent (put a letter ``D`` in the column ``FCD`` of the FTBL sections ``FLUXES/NET`` and ``FLUXES/XCH``) and see an error message that will suggest candidates for free fluxes. Another option is to use ``--ffguess`` flag that will automatically partition not constrained fluxes between free and dependent.

 * if your resulting fluxes are badly defined (statistically or structurally), i.e. they have big confidence intervals or the Jacobian is rank deficient, you can try to play with input labeling (cf. IsoDesign software at http://metatoul.insa-toulouse.fr/metasys/software/isodes/) or try to collect some additional data on metabolites not yet measured. To have some insights on what part of the network is already well defined and which one still needs additional measurements, you can try to run ``influx_si`` with an option ``--ln`` (as `least norm`) (in addition to ``--noopt`` option) and examine standard deviation of the fluxes/concentrations in the resulting KVH file. Another possibility is to use `parallel labeling experiments <https://doi.org/10.1016/j.ymben.2012.11.010>`_ (cf. manual section :ref:`Parallel experiments <prlexp>`)

Once again, if you could not resolve your problem during these steps, see the next section  :doc:`Consulting and more <consulting>`.
