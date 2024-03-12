# xtb_ase_io_calculator
`ASE` `FileIOCalculator` child class for the `xTB` code of Grimme and coworkers. 

The code in this repository was used, coupled with the GFN0-xTB Hamiltonian in the work described in: 
- Veronika Jurásková, Frederic Célerse, Ruben Laplaza, Clemence Corminboeuf; Assessing the persistence of chalcogen bonds in solution with neural network potentials. J. Chem. Phys. 21 April 2022; 156 (15): 154112. https://doi.org/10.1063/5.0085153

`ASE` (Atomic Simulation Environment) can be found in:
- https://wiki.fysik.dtu.dk/ase/
- https://gitlab.com/ase/ase

The underlying `xTB` code can be found in:
- https://xtb-docs.readthedocs.io
- https://github.com/grimme-lab/xtb

The calculator requires `ASE`, a working `xTB` install and the environment variable `$XTBHOME` to be set.
It also has explicit dependencies on `numpy` and `subprocess`.
