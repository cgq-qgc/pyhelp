copy %PREFIX%\libs\libpython%CONDA_PY%.dll.a %PREFIX%\libs\libpython%CONDA_PY%.a
"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt
if errorlevel 1 exit 1
