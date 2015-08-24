#!/usr/bin/env python

"""
The script create a setup folder to run Yank with a articular receptor and ligand.

The script uses a simple heuristic algorithm to ensure that the receptor and the
ligand are very close and do not overlap at the beginning of the simulation.

Examples
--------
Setup Yank to run Abl in the configuration described by 3UE4 and bosutinib in
implicit solvent.

> python setup_yank.py --receptor abl/3UE4-pdbfixer.pdb --ligand Bosutinib.mol2

Setup Yank to run Abl in the configuration described by 2HYY and bosutinib in
explicit solvent.

> python setup_yank.py --receptor abl/2HYY-pdbfixer.pdb --ligand Imatinib.mol2 --solvate

"""


import os
import shutil
import subprocess
from contextlib import contextmanager

import numpy as np
from openeye import oechem


# Constants and configuration
#----------------------------
SETUP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'setup')
KINASES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kinases')
LIGANDS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ligands')
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils')

LEAP_IN_IMPLICIT = os.path.join(UTILS_DIR, 'setup.leap-implicit.in')
LEAP_IN_EXPLICIT = os.path.join(UTILS_DIR, 'setup.leap-explicit.in')

# Min and max distance of the ligand at the beginning of the simulation
MIN_DISTANCE = 1.5
MAX_DISTANCE = 2.5


# Utility functions
#----------------------------
@contextmanager
def working_directory(path):
    """Context to set the working directory to the given path."""
    current_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(current_dir)


def run_command(cmd):
    """Run a command with subprocess and raise RuntimeError if execution fails."""
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        raise RuntimeError('{!r} failed, exit status: {:d}'.format(cmd, return_code))


def pull_close(receptor_oemol, ligand_oemol, min_bound, max_bound):
    """Heuristic algorithm to quickly translate the ligand close to the receptor.

    The distance of the ligand from the receptor here is defined as the shortest
    Euclidean distance between an atom of the ligand and one of the receptor.

    The molecules positions will not be modified if the ligand is already at a
    distance in the interval [min_bound, max_bound].

    Parameters
    ----------
    receptor_oemol : openeye.oechem
        The receptor molecule.
    ligand_oemol : openeye.oechem
        The ligand molecule to push close to the receptor.
    min_bound : float
        Minimum distance from the receptor to the ligand. This should be high
        enough for the ligand to not overlap the receptor atoms at the beginning
        of the simulation.
    max_bound : float
        Maximum distance from the receptor to the ligand. This should be short
        enough to make the ligand and the receptor interact since the beginning
        of the simulation.

    """

    goal_distance = (min_bound + max_bound) / 2

    # Get all the ligand and receptor atoms position
    oe_coords = oechem.OEFloatArray(3)

    ligand_pos = np.zeros((ligand_oemol.NumAtoms(), 3))
    for i, atom in enumerate(ligand_oemol.GetAtoms()):
        ligand_oemol.GetCoords(atom, oe_coords)
        ligand_pos[i] = oe_coords

    receptor_pos = np.zeros((receptor_oemol.NumAtoms(), 3))
    for i, atom in enumerate(receptor_oemol.GetAtoms()):
        receptor_oemol.GetCoords(atom, oe_coords)
        receptor_pos[i] = oe_coords

    # Find translation
    final_translation = np.zeros(3)
    while True:

        # Compute squared distances
        # Each row is an array of distances from a ligand atom to all receptor atoms
        # We don't need to apply square root to everything
        distances2 = np.array([((receptor_pos - ligand_pos[i])**2).sum(1)
                               for i in xrange(len(ligand_pos))])

        # Find closest atoms and their distance
        min_idx = np.unravel_index(distances2.argmin(), distances2.shape)
        min_dist = np.sqrt(distances2[min_idx])

        # If closest atom is between boundaries translate ligand
        if min_bound <= min_dist <= max_bound:
            break

        # Compute unit vector that connects receptor and ligand atom
        direction = receptor_pos[min_idx[1]] - ligand_pos[min_idx[0]]
        direction = direction / np.sqrt((direction**2).sum())

        if max_bound < min_dist:  # the atom is far away
            translation = (min_dist - goal_distance) * direction
            ligand_pos += translation
            final_translation += translation
        elif min_dist < min_bound:  # the two molecules overlap
            max_dist = np.sqrt(distances2.max())
            translation = (max_dist + goal_distance) * direction
            ligand_pos += translation
            final_translation += translation

    # Translate OEChem molecule
    translation_oe = oechem.OEDoubleArray(final_translation)
    oechem.OETranslate(ligand_oemol, translation_oe)


# Main function
#----------------------------
def main(receptor_file_name, ligand_file_name, setup_directory_name, solvate=False):

    # Cleanup setup directory
    if os.path.exists(setup_directory_name):
        shutil.rmtree(setup_directory_name)
    os.makedirs(setup_directory_name)

    # Read ligand and receptor molecule
    ifs_mol2 = oechem.oemolistream()
    ifs_mol2.open(ligand_file_name)
    ligand_oemol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs_mol2, ligand_oemol)
    ifs_mol2.close()

    ifs_mol2 = oechem.oemolistream()
    ifs_mol2.open(receptor_file_name)
    receptor_oemol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs_mol2, receptor_oemol)
    ifs_mol2.close()

    # Push ligand close to receptor
    pull_close(receptor_oemol, ligand_oemol, MIN_DISTANCE, MAX_DISTANCE)

    # Add residue name 'MOL'
    residue = oechem.OEResidue()
    residue.SetName('MOL')
    for atom in ligand_oemol.GetAtoms():
        oechem.OEAtomSetResidue(atom, residue)

    # Parametrize ligand
    with working_directory(setup_directory_name):

        # Save translated ligand
        ofs = oechem.oemolostream()
        ofs.open('ligand.mol2')
        oechem.OEWriteMolecule(ofs, ligand_oemol)
        ofs.close()

        # Parametrize ligand
        print "Parameterizing ligand with GAFF..."
        run_command('antechamber -fi mol2 -i ligand.mol2 -fo mol2 -o ligand.gaff.mol2')
        run_command('parmchk -i ligand.gaff.mol2 -o ligand.gaff.frcmod -f mol2')

        # Copy receptor so that leap will know the PDB file name
        shutil.copyfile(receptor_file_name, 'receptor.pdb')

        # Create AMBER prmtop/inpcrd files.
        print "Creating AMBER prmtop/inpcrd files..."
        cmd = 'tleap -f {!s} > setup.leap.out'
        if solvate:
            run_command(cmd.format(LEAP_IN_EXPLICIT))
        else:
            run_command(cmd.format(LEAP_IN_IMPLICIT))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--receptor', metavar='RECEPTORPDB', required=True,
                        help='The name of the receptor PDB file')
    parser.add_argument('-l', '--ligand', metavar='LIGANDMOL2', required=True,
                        help='The name of the charged ligand MOL2 file')
    parser.add_argument('-s', '--solvate', action='store_true')
    args = parser.parse_args()

    receptor_file_name = os.path.join(KINASES_DIR, args.receptor)
    ligand_file_name = os.path.join(LIGANDS_DIR, args.ligand)
    setup_directory_name = os.path.join(SETUP_DIR, args.receptor.split('.')[0], args.ligand.split('.')[0])

    main(receptor_file_name, ligand_file_name, setup_directory_name, solvate=args.solvate)
