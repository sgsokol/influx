if (!require("devtools")) install.packages("devtools", repos="https://cloud.r-project.org/")
devtools::update_packages(c("nnls", "rmumps", "arrApply", "slam", "limSolve", "devtools"))
devtools::install_github("sgsokol/multbxxc")
