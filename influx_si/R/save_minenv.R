# save minimal environement for RPP simulation
if (write_res)
    save(param, labargs, file=file.path(dirres, "tmp", sprintf("%s.RData", baseshort)))
