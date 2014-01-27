
.. _howto:

==========
How to ...
==========

.. describe:: ... choose free fluxes?

 You can define in FTBL all not constrained fluxes as dependent (put a letter ``D`` in the column ``FCD`` of the FTBL sections ``FLUXES/NET`` and ``FLUXES/XCH``), run ``influx_s`` and see an error message that will suggest some candidates for free fluxes. For these fluxes, put a letter ``F`` in the column ``FCD`` and some numeric value in the next column ``VALUE(F/C)`` to provide a starting value for the fitting. Don't use 0 as starting value as it might lead to singular matrices in cumomer balances.

.. describe:: ... get statistical information for a given set of free fluxes without
   fitting measurements?

 Put these values in the corresponding FTBL file as starting values for free fluxes and use ``influx_s`` with ``--noopt`` option.

.. describe:: ... accelerate calculations?

 You can relax stopping criterion and pass from 1.e-5 (by default) to, for example, 1.e-2 if this precision is sufficient for you. Use ``optctrl_errx`` option in FTBL file (section ``OPTIONS``) for this.

 If you mean to accelerate Monte-Carlo simulations in Unix environment, you can use a hardware with many cores. In this case, the wall clock time can be reduced significantly. Note that distant nodes, even inside of the same cluster, are not used in the such kind of Monte-Carlo simulations.

 Check that your system is not using swap (disk) memory. If it is the case, stop other applications running in parallel with ``influx_s``. If possible extend the RAM on your hardware.

.. describe:: ... extend upper limit for non linear iterations?

 By default, this value is 50 which should be largely sufficient for most cases. If not, you can set another value via ``optctrl_maxit`` option in the FTBL file (section ``OPTIONS``). But most probably, you would like to check your network definition or to add some data or to change a substrate labeling, anyway to do something to get a well defined network instead of trying to make converge the fitting on some biologically almost meaningless situation.