
.. _trouble:

===============
Troubleshooting
===============

The software is provided "AS IS" so for the troubleshooting you are on your own. We don't provide any support of any kind.

Anyway, if you have a problem during installation, you can ask for help from your local computer desk.

If you have a problem with FTBL editing, you can read the documentation from `13CFlux <https://www.13cflux.net>`_ and/or interpret error messages generated during FTBL parsing.

If you have some difficulties in choosing free fluxes, define all not constrained fluxes as dependent (put a letter ``D`` in the column ``FCD`` of the FTBL sections ``FLUXES/NET`` and ``FLUXES/XCH``) and see an error message that will suggest candidates for free fluxes.

Finally, if your resulting fluxes are badly statistically or structurally defined, i.e. they have big confidence intervals or the Jacobian is rank deficient, you can try to play with input labeling or try to collect some additional data on metabolites not yet measured. To have some insights on what part of the network is already well defined and which one still need additional measurements, you can try to run influx_s with an option ``--ln`` (as `least norm`) and examine standard deviation of the fluxes in the resulting KVH file.
