# 2021-02-25 sokol@insa-toulouse.fr

FTBL files for single labeling experiment
-----------------------------------------
 * e_coli.ftbl -- real-world example on central carbon metabolism of E.coli (authored by P. Millard)
    run with: influx_s e_coli
 * e_coli_growth.ftbl -- same as e_coli.ftbl but with growth fluxes modeled as -µ*[metab_concentration]
    run with: influx_s e_coli_growth
 * e_coli_i.ftbl -- instationary labeling in E.coli on synthetic data.
    run with: influx_i e_coli_i
 * e_coli_iv.ftbl -- instationary labeling in E.coli with input labeling in a form of periodic rectangular pulses (no fitting, only simulation)
    run with: influx_i e_coli_iv
 * ex_i_2box_var.ftbl -- instationary labeling in a linear metabolic chain of 2 metabolites with label input in an exponential form (no fitting, only simulation)
    run with: influx_i ex_i_2box_var

FTBL files for parallel labeling experiments
---------------------------------------------
 * prl_exp/e_coli_glc1-6n.ftbl -- 6 parallel experiments, each having as label input a Glucose molecule labeled in corresponding atom. 
    run with: influx_s prl_exp/e_coli_glc1-6n
 * prl_exp/e_coli_GX_prl.ftbl -- 2 parallel experiments with instationary labeling
    run with: influx_i prl_exp/e_coli_GX_prl

Helper files
------------
 * README.txt -- this file
 * e_coli_iv_funlab.R - helper for e_coli_iv.ftbl, defines some variables referenced in LABEL_INPUT R expressions
 * e_coli_msen.txt -- kinetic labeling data for e_coli_i.ftbl
 * prl_exp/e_coli_glc2n.ftbl, prl_exp/e_coli_glc3n.ftbl, prl_exp/e_coli_glc4n.ftbl, prl_exp/e_coli_glc5n.ftbl, prl_exp/e_coli_glc6n.ftbl -- auxiliary FTBL files for parallel experiments referenced from prl_exp/e_coli_glc1-6n.ftbl
 * prl_exp/e_coli_GX_X.ftbl -- auxiliary FTBL for parallel experiment referenced from prl_exp/e_coli_GX_prl.ftbl
 * prl_exp/e_coli_GX_MS.txt, prl_exp/e_coli_GX_X_MS.txt -- simulated kinetic labeling data for parallel experiments referenced from prl_exp/e_coli_GX_prl.ftbl and prl_exp/e_coli_GX_X.ftbl
 * ok/* -- results of all cited FTBL files, to be used as reference in comparison with user's results
