@echo on
rem echo %%0 = %0
set direx=%~dp0
set DEBUG=
if (%2) == (DEBUG) set DEBUG=DEBUG
set me=%~nx0

set org=%~dpn1
if not "dummy%org%" == "dummy" goto :work
echo %me%: expecting ftbl file name
set ERR=1
goto :end

:work
python "%direx%influx_i.py" "%org%" %DEBUG% >"%org%.err"
R CMD SHLIB --preclean "%org%.f"
R --no-save --silent --args --meth %DEBUG% < "%org%.R" > "%org%.log"

:end
rem the end
