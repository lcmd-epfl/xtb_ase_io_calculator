# xtb_ase_io_calculator
'ASE' 'FileIOCalculator' child class for the 'xTB' code of Grimme and coworkers. 

The code in this repository was used, coupled with the GFN0-xTB Hamiltonian in the work described in: 
- chemrXiv preprint 

ASE (Atomic Simulation Environment) can be found in:
- https://wiki.fysik.dtu.dk/ase/
- https://gitlab.com/ase/ase

The underlying xTB code can be found in:
- https://xtb-docs.readthedocs.io
- https://github.com/grimme-lab/xtb

The calculator requires 'ASE', a working 'xTB' install and the environment variable '$XTBHOME' to be set.
It also has explicit dependencies on 'numpy' and 'subprocess'.
