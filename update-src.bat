@echo off
echo ---------------------------------------------------
echo Dr. Probe command-line tool - Source file update.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
SET LPATH=F:\VirtualShare\drprobe_clt\
SET FPATH=F:\Data\F90\drprobe_clt\
echo Updating source files in [%LPATH%]
echo from base directory [%FPATH%].
copy "%FPATH%*.txt" "%LPATH%" /Y
copy "%FPATH%*.md" "%LPATH%" /Y
copy "%FPATH%*.sh" "%LPATH%" /Y
copy "%FPATH%LICENSE" "%LPATH%" /Y
echo - subfolder [common] ...
copy "%FPATH%common\*.txt" "%LPATH%common\" /Y
copy "%FPATH%common\*.f" "%LPATH%common\" /Y
copy "%FPATH%common\*.in" "%LPATH%common\" /Y
copy "%FPATH%common\*.f90" "%LPATH%common\" /Y
echo - subfolder [celslc] ...
copy "%FPATH%celslc\*.txt" "%LPATH%celslc\" /Y
copy "%FPATH%celslc\*.f" "%LPATH%celslc\" /Y
copy "%FPATH%celslc\*.f90" "%LPATH%celslc\" /Y
copy "%FPATH%celslc\makefile" "%LPATH%celslc\" /Y
echo - subfolder [msa] ...
copy "%FPATH%msa\*.txt" "%LPATH%msa\" /Y
copy "%FPATH%msa\*.f" "%LPATH%msa\" /Y
copy "%FPATH%msa\*.f90" "%LPATH%msa\" /Y
copy "%FPATH%msa\makefile" "%LPATH%msa\" /Y
echo - subfolder [wavimg] ...
copy "%FPATH%wavimg\*.txt" "%LPATH%wavimg\" /Y
copy "%FPATH%wavimg\*.f" "%LPATH%wavimg\" /Y
copy "%FPATH%wavimg\*.f90" "%LPATH%wavimg\" /Y
copy "%FPATH%wavimg\makefile" "%LPATH%wavimg\" /Y
echo ---------------------------------------------------
echo Done.