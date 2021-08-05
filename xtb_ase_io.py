#!/usr/bin/env python
"""This module defines an ASE interface to XTB 6.2.0 through I/O."""

import os
import numpy as np
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.constraints import FixAtoms
from ase.units import Bohr, Hartree
from subprocess import check_output, STDOUT

# $XTBPATH must point to the XTB installation directory


def get_xtb_version(command, verb=0):
    """Executes command using subprocess and tries to find the xTB version from the STDOUT."""
    check = False
    for line in (
        check_output(command, shell=True, stderr=STDOUT).decode("ascii").split("\n")
    ):
        if "VERSION" in line.upper():
            check = True
            for i, word in enumerate(line.upper().split(" ")):
                if "VERSION" in word:
                    version = line.split(" ")[i + 1]
    if not check:
        raise RuntimeError("Cannot recognize XTB version")
    if verb > 0:
        print(f"XTB version {version}")
    return version


class XTB(FileIOCalculator):
    """A simple FileIOCalculator child class that implements energies and forces using I/O from xTB.
    Importantly, requires that the environment variable $XTBPATH is set and that command is the effective
    command line execution argument of xTB. Charge and spin are inherited from the Atoms object that the
    calculator will be atteched to in ASE.

    Parameters:
    -----------
    :param restart: Passes the restart flag to the FileIOCalculator base class.
    :param label: Label that will be used to create a directory to run calculations.
    :param atoms: An ASE Atoms object that defines charge and multiplicity.
    :param command: A string corresponding to the xTB execution command.
    :param verb: A verbosity flag that controls printing to STDOUT from this calculator.
    """

    implemented_properties = ["energy", "forces"]
    discard_results_on_any_change = True

    def __init__(
        self,
        restart=None,
        label="xtbcalc/",
        atoms=None,
        command="xtb --gfn 0 --grad  xtb.tmol > xtb.log",
        verb=0,
    ):
        xtbpath = os.environ.get("XTBPATH")
        if xtbpath is None:
            mess = "The XTBPATH variable is not defined."
            raise ValueError(mess)
        self.label = label
        self.verb = verb
        FileIOCalculator.__init__(
            self,
            restart=restart,
            label=label,
            atoms=atoms,
            command=command,
        )
        self.clean()
        get_xtb_version(command + " --dry --version 1>&2 ", self.verb)

    def write_charge_mult(self):
        charge = int(self.atoms.get_initial_charges().sum().round())
        uhf = int(self.atoms.get_initial_magnetic_moments().sum().round())
        with open(self.label + ".CHRG", "w") as f1:
            f1.write(str(charge))
        with open(self.label + ".UHF", "w") as f2:
            f2.write(str(uhf))

    def write_input(self, atoms, properties=None, system_changes=None):
        if not (np.all(atoms.pbc) or not np.any(atoms.pbc)):
            raise RuntimeError("PBC must be all true or all false")
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        filename = self.label + "xtb.tmol"
        coord = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        self.write_charge_mult()
        with open(filename, "w") as fd:
            fix_indices = set()
            if self.atoms.constraints:
                for constr in atoms.constraints:
                    if isinstance(constr, FixAtoms):
                        fix_indices.update(constr.get_indices())

            fix_str = []
            for i in range(len(atoms)):
                if i in fix_indices:
                    fix_str.append("f")
                else:
                    fix_str.append("")
            if np.all(atoms.pbc):
                fd.write("$periodic 3\n")
                fd.write("$cell angs\n")
                fd.write(
                    " {0} {1} {2}  ".format(
                        atoms.get_cell_lengths_and_angles()[0],
                        atoms.get_cell_lengths_and_angles()[1],
                        atoms.get_cell_lengths_and_angles()[2],
                    )
                )
                fd.write(
                    " {0} {1} {2}\n".format(
                        atoms.get_cell_lengths_and_angles()[3],
                        atoms.get_cell_lengths_and_angles()[4],
                        atoms.get_cell_lengths_and_angles()[5],
                    )
                )
            fd.write("$coord angs\n")
            for (x, y, z), s, fix in zip(coord, symbols, fix_str):
                fd.write(
                    "%20.14f  %20.14f  %20.14f      %2s  %2s \n"
                    % (x, y, z, s.lower(), fix)
                )
            fd.write("$end\n")

    def read_results(self):
        self.read_energy(self.directory, self.label, self.verb)
        self.read_forces(self.atoms, self.directory, self.label, self.verb)

    def read_forces(self, atoms, directory, label, verb=0):
        gradient_file = label + "gradient"
        natoms = len(atoms)
        if not os.path.exists(gradient_file):
            raise ReadError("The gradient file does not exist.")
        if os.path.isfile(gradient_file):
            with open(gradient_file, "r") as fd:
                lines = fd.readlines()
        lines = lines[-1 - natoms : -1]
        assert len(lines) == natoms
        self.results["forces"] = np.zeros((natoms, 3), float)
        for i in range(natoms):
            line = [s for s in lines[i].strip().split(" ") if len(s) > 0]
            f = -np.array([float(x.replace("D", "E")) for x in line[0:3]])
            self.results["forces"][i, :] = f * (Hartree / Bohr)
        os.remove(gradient_file)
        if verb > 2:
            print("Forces of the last atom:\n {0}".format(f * (Hartree / Bohr)))

    def read_energy(self, directory, label, verb=0):
        energy_file = label + "energy"
        if not os.path.exists(energy_file):
            raise ReadError("The energy file does not exist.")
        if os.path.isfile(energy_file):
            with open(energy_file, "r") as fd:
                lines = fd.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("$end"):
                escf = (
                    float(lines[i - 1].strip().split(" ")[4].replace("D", "E"))
                    * Hartree
                )
                self.results["energy"] = escf
                os.remove(energy_file)
                break
        if verb > 1:
            print("Energy:\n {0}".format(escf))

    def finished_successfully(self):
        finished = False
        mess = ""
        for line in open(self.label + "xtb.log", "r"):
            if line.rfind("finished run") > -1:
                finished = True
            if line.contains("error"):
                mess += line
        if self.verb > 0:
            if finished:
                print("Normal termination.")
            else:
                print("Error termination.")
        return finished, mess

    def clean(self):
        files_to_clean = [
            self.label + "xtb.log",
            self.label + "energy",
            self.label + "gradient",
            self.label + "xtbrestart",
            self.label + "charges",
        ]
        for f in files_to_clean:
            try:
                os.remove(f)
            except OSError:
                pass


def test_xtb_ase_io():
    from ase.build import molecule

    atoms = molecule("O2")
    atoms.set_calculator(XTB(command="xtb --gfn 2 --grad  xtb.tmol > xtb.log", verb=1))
    print("Potential energy %5.5f eV" % atoms.get_potential_energy())
    print(
        "Forces per atom:\n {0} \n {1}".format(
            atoms.get_forces()[0][:], atoms.get_forces()[1][:]
        )
    )


if __name__ == "__main__":
    test_xtb_ase_io()
