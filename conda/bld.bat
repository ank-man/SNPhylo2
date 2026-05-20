@echo off
:: Windows build script for SNPhylo2

%PYTHON% -m pip install . --no-deps --ignore-installed -vv
if errorlevel 1 exit 1

echo SNPhylo2 installation complete
