
.. _trouble:

===============
Troubleshooting
===============

The software is provided "AS IS" so for the troubleshooting you are on your own. We don't provide any support of any kind.

Anyway, if you have a problem during installation, you can ask for help from your local computer desk.

If you have a problem with FTBL editing, you can read the documentation from `13CFlux <https://www.13cflux.net>`_ and/or interpret error messages generated during FTBL parsing.

If you have some difficulties in choosing free fluxes, define all not constrained fluxes as dependent (put a letter ``D`` in the column ``FCD`` of the FTBL sections ``FLUXES/NET`` and ``FLUXES/XCH``) and see an error message that will suggest candidates for free fluxes.

If your resulting fluxes are badly statistically or structurally defined, i.e. they have big confidence intervals or the Jacobian is rank deficient, you can try to play with input labeling (cf. IsoDesign software at http://metasys.insa-toulouse.fr/software/isodes/) or try to collect some additional data on metabolites not yet measured. To have some insights on what part of the network is already well defined and which one still needs additional measurements, you can try to run influx_s with an option ``--ln`` (as `least norm`) (in addition to ``--noopt`` option) and examine standard deviation of the fluxes in the resulting KVH file.

Finaly, if you think to discover a bug in ``influx_s`` you can report it to the author by email to ``sokol [at] insa-toulouse [dot] fr``. At this moment, please be sure to use the latest available release as the bug may be already corrected or not actual any more. Note also that we can't guarantee that any particular bug can be fixed in any particular release or can be fixed at all.