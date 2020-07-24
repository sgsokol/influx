
.. _howto:

==========
How to ...
==========

.. describe:: ... choose free fluxes?

 You can define in FTBL all not constrained fluxes as dependent (put a letter ``D`` in the column ``FCD`` of the FTBL sections ``FLUXES/NET`` and ``FLUXES/XCH``), run ``influx_si`` and see an error message that will suggest some candidates for free fluxes. For these fluxes, put a letter ``F`` in the column ``FCD`` and some numeric value in the next column ``VALUE(F/C)`` to provide a starting value for the fitting. Don't use 0 as starting value as it might lead to singular matrices in cumomer balances.
 
 If you want to create an FTBL *de novo*, consider using application ``txt2ftbl.py`` included in ``influx_si`` package. Not only it translates an easily readable/writable text format into FTBL one, but it also automatically assigns some fluxes to be free.

.. describe:: ... get statistical information for a given set of free fluxes without
   fitting measurements?

 Put these values in the corresponding FTBL file as starting values for free fluxes and use ``influx_si`` with ``--noopt`` option.

.. describe:: ... accelerate calculations?

 You can relax stopping criterion and pass from 1.e-5 (by default) to, for example, 1.e-2 if this precision is sufficient for you. Use ``optctrl:nlsic:errx`` option in FTBL file (section ``OPTIONS``) for this.

 If you mean to accelerate Monte-Carlo simulations in Unix environment, you can use a hardware with many cores. In this case, the wall clock time can be reduced significantly. Note that distant nodes, even inside of the same cluster, are not used in the such kind of Monte-Carlo simulations.

 Check that your system is not using swap (disk) memory. If it is the case, stop other applications running in parallel with ``influx_si``. If possible extend the RAM on your hardware.

.. describe:: ... extend upper limit for non linear iterations?

 By default, this value is 50 which should be largely sufficient for most cases. If not, you can set another value via ``optctrl:nlsic:maxit`` option in the FTBL file (section ``OPTIONS``). But most probably, you would like to check your network definition or to add some data or to change a substrate labeling, anyway to do something to get a well defined network instead of trying to make converge the fitting on some biologically almost meaningless situation.

.. describe:: ... make FTBL file with synthetic data?

 Follow for example steps outlined hereafter:
  - edit FTBL file(s) with ``NA`` in measurements and realistic SD, name it e.g. ``new_NA.ftbl``
  - simulate data: ::
  
    $ influx_s.py --noopt --addnoise new_NA
    
  - prepare FTBL sections with simulated data: ::
  
     $ res2ftbl_meas.py new_NA_res.kvh
    
    It will create file (or files if there are parallel experiments) with synthetic data formatted for inclusion in FTBL file: ``new_NA_sim1.ftbl``, ``new_NA_sim2.ftbl``, etc.)
  - copy/paste simulated data to a new file ``new.ftbl`` with numeric data instead of ``NA``.
  - use FTBL with synthetic data: ::
  
    $ influx_s.py new.ftbl

.. describe:: ... do custom post-treatment of Monte-Carlo iterations?

 Let suppose you want to filter some of Monte-Carlo (MC) iterations based on their cost values.
 In ``OPTIONS/posttreat_R`` of your FTBL file add ``save_all.R``. The file ``save_all.R`` can be found in ``R`` directory of ``influx_si`` distribution. Execution of ``save_all.R`` at the end of calculations will simply save all session variables in ``mynetwork.RData`` file (supposing that your FTBL file is named ``mynetwork.ftbl``). In particular, you need ``free_mc`` matrix which contains free parameters (each column results from a given MC iteration). After that you can open an interactive R session in your working directory and run something similar to:
 
 .. code-block:: r
  
  # preparations
  load("mynetwork.RData")
  source(file.path(dirx, "libs.R"))
  source(file.path(dirx, "opt_cumo_tools.R"))
  #source(file.path(dirx, "opt_icumo_tools.R")) # uncoment for influx_i use
  tmp=sparse2spa(spa)
  
  # doing something useful
  # here, we calculate a vector of cost values, one per MC iteration
  cost_mc=apply(free_mc, 2, function(p) cumo_cost(p, labargs))
  # do something else ...

 If, instead of cost values, you need for example a full set of net-xch fluxes then do
 
 .. code-block:: r
  
  allflux_mc=apply(free_mc, 2, function(p) param2fl(p, labargs)$fallnx)
   
 for residuals, do:
 
 .. code-block:: r
  
  resid_mc=apply(free_mc, 2, function(p) lab_resid(p, FALSE, labargs)$res)

 After that, you can filter or do whatever needed with obtained vectors and matrices.


