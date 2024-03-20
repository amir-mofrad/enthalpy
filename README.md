# Enthalpy of Formation--VASP
The main purpose of this Python code is to calculate the enthalpy of formation from DFT calculations that are done using the VASP code.
Other structural parameters are also calculated (i.e., cell parameters, coordination number, etc.). However, the sole purpose was to
have a code that quickly (not really sure about that) and easily calculates the enthalpy of formation (from elements) of a compound.
A few things to consider when using this code:
  1. The fitted atomic potentials are taken from the Open Quantum Materials Database (OQMD) paper, fit-all column. https://doi.org/10.1038/npjcompumats.2015.10
  2. When using/calling this code make sure you have the OUTCAR and either POSCAR or CONTCAR.
  3. Change/remove the first line depending on your Python environment.
  4. Your planewave cutoff energy should be 520 eV to be compatible with OQMD.
  5. The atomic potentials of actinides were calculated from our previous work.
  6. Feel free to uncomment print commands for your needs, i.e., kJ/mol vs eV/atom. 
  7. Also, yes, this is 0 K enthalpy of formation.
  
Good Luck!
