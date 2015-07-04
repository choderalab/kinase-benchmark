#!/usr/bin/env python

"""
This script reads all the ligands in SMILES format contained in the CSV file
indicated by the constant SMILES_FILE_PATH. Then it uses the OpenEye Omega
toolkit to assign a structure to the molecule and Schrodinger epik to assign
its protonation state. The script assigns canonical AM1-BCC charges to the
ligand as well and saves the final molecule in the folder LIGANDS_DIR_PATH.

It requires OpenEye toolkit and Schrodinger installed. Moreover the environment
variable $SCHRODINGER must point to its installation directory (on Mac OS X it
is usually '/opt/schrodinger/suites201X').
"""

import os
import shutil
import tempfile
from contextlib import contextmanager
from openeye import oechem, oeomega
from openmoltools.openeye import smiles_to_oemol
from utils.tpl_maker_am1bcc import mk_single_conformer_epik


# Constants and configuration
#----------------------------
LIGANDS_DIR_PATH = os.path.join(os.path.dirname(__file__), 'ligands')
SMILES_FILE_PATH = os.path.join(os.path.dirname(__file__), 'clinical-kinase-inhibitors.csv')

PH = 7  # pH to pass to epik to generate protonation states
OMEGA_MAX_CONFS = 1  # How many conformers Omega must generate


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


if __name__ == '__main__':

    # Make directory for temporary epik output
    temp_dir = tempfile.mkdtemp()

    # Create output directory
    if not os.path.exists(LIGANDS_DIR_PATH):
        os.makedirs(LIGANDS_DIR_PATH)

    # Parse SMILES file and generate ligand
    smiles_file = open(SMILES_FILE_PATH, 'r')
    for line in smiles_file:

        # Get SMILES representation from CVS file
        ligand_name, smiles_str = line.strip().split(',')

        print "Generating inhibitor", ligand_name

        # Convert to OEMol molecule
        mol_oemol = smiles_to_oemol(smiles_str)
        mol_oemol.SetTitle(ligand_name)

        # Use OpenEye Omega toolkit to generate lowest energy structure
        omega = oeomega.OEOmega()
        omega.SetCanonOrder(False)
        omega.SetMaxConfs(OMEGA_MAX_CONFS)
        omega(mol_oemol)

        # Generate protonation/tautomeric states with epik and charge molecule
        valid_structure = False
        with working_directory(temp_dir):
            for i in range(mol_oemol.GetMaxConfIdx()):
                try:
                    confomer = mol_oemol.GetConf(oechem.OEHasConfIdx(i))
                    mk_single_conformer_epik(ligand_name, confomer, pH=PH)
                    valid_structure = True
                    break
                except ValueError:
                    print "WARNING: Omega structure #{:d} discarded.".format(i)

        # Copy to output directory final file
        if valid_structure:
            mol2_file_name = ligand_name + '.mol2'
            src_path = os.path.join(temp_dir, mol2_file_name)
            dst_path = os.path.join(LIGANDS_DIR_PATH, mol2_file_name)
            shutil.copyfile(src_path, dst_path)

    smiles_file.close()
