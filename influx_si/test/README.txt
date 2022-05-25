# 2022-05-25 sokol@insa-toulouse.fr

MTF files for single labeling experiment
----------------------------------------
 * mtf/e_coli -- real-world example on central carbon metabolism of E.coli (original e_coli.ftbl authored by P. Millard)
    run with: influx_s --prefix e_coli
 * mtf/e_coli_growth.* -- same as e_coli but with growth fluxes modeled as -µ*[specie_concentration]
    run with: influx_s --prefix e_coli_growth
 * mtf/e_coli_i -- instationary labeling in E.coli on synthetic data.
    run with: influx_i --prefix e_coli_i
 * mtf/e_coli_iv -- instationary labeling in E.coli with input labeling in a form of periodic rectangular pulses (no fitting, only simulation)
    run with: influx_i --prefix e_coli_iv
 * mtf/ex_i_2box_var -- instationary labeling in a linear metabolic chain of 2 metabolites with label input in an exponential form (no fitting, only simulation)
    run with: influx_i --prefix ex_i_2box_var

MTF files for parallel labeling experiments
-------------------------------------------
 * prl_exp/mtf/e_coli_glc1-6n.ftbl -- 6 parallel experiments, each having as label input a Glucose molecule labeled in corresponding atom. 
    run with: influx_s --prefix prl_exp/mtf/e_coli_glc1-6n
 * prl_exp/mtf/e_coli_GX_prl -- 2 parallel experiments with instationary labeling
    run with: influx_i --prefix prl_exp/e_coli_GX_prl

Helper files
------------
 * README.txt -- this file
 * mtf/e_coli_iv_funlab.R - helper for e_coli_iv.*, defines some variables referenced in .linp R expressions
 * prl_exp/mtf/e_coli_glc2n.*, prl_exp/mtf/e_coli_glc3n.*, prl_exp/mtf/e_coli_glc4n.*, prl_exp/mtf/e_coli_glc5n.*, prl_exp/mtf/e_coli_glc6n.* -- auxiliary .miso nad .linp files for parallel experiments referenced from prl_exp/mtf/e_coli_glc1-6n.opt
 * prl_exp/mtf/e_coli_GX_X.* -- auxiliary .miso, .linp and .opt files for parallel experiment referenced from prl_exp/mtf/e_coli_GX_prl.opt
 * ok/* -- results of all cited MTF files, to be used as reference in comparison with user's results
