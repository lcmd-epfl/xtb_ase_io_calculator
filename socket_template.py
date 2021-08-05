#!/usr/bin/env python3

import numpy as np
from xtb_ase_io import XTB
from ase.io import read
from ase.calculators.socketio import SocketClient

path = "example.xyz"
host = "127.127.127.127"
port = "0000"
mol = read(path, format="xyz")

atnums = mol.get_atomic_numbers()
coords = mol.get_positions()

# Charge must be set manually as with any ASE method, only the global charge is taken into account for XTB
charges = np.zeros(len(atnums))
charges[0] = -1
mol.set_initial_charges(charges)
cell = mol.get_cell()
pbc = mol.get_pbc()

# Calculator setup
mol.set_calculator(XTB())
client = SocketClient(host=host, port=port)
client.run(mol)
