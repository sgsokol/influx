# This file can be used as a postreatment script in influx_s
# software (cf. parameter posttreat_R in FTBL/OPTIONS section)

# Purpose: save the whole environement in a mynetwork.RData file
# to explore later in an interactive R session.

# file name to save
f=file.path(dirw, sprintf("%s.RData", baseshort))
save(list = ls(all = TRUE), file = f)
