#!/usr/bin/python

"""
tpl_maker.py created by Denise M Kilburg for the Gunner Lab at CUNY-City College for the
generation of ligand topology files for use in MCCE.

tpl_maker.py requires a protein data bank file, ligand code, chain id, and charge.

for example: tpl_maker.py 4Q6T.pdb PO4 A -3

At the moment, the charge represents the most negative charge on an ion and it dictates
the amount of conformers made. For example: in the above PO4 example, charge = -3 will
create a PO401(neutral conformer), PO4-1, PO4-2,PO4-3 conformers. At the moment, the
charge also specifies the number of hydrogens to be removed from the neutral conformer.
At present, the program removes hydrogens from hydroxyl groups. If there are no
hydroxyl groups it will remove hydrogens starting from H1, H2,...
The add_hydrogens function adds hydrogens in 2 ways. If an atom is connected to 2 or more
atoms, it looks at the angles that subtend them and assigns orbitals based on geometry:
sp3 == tetrahedral, sp2 == trigonal planar, sp == linear. It then adds the proper number of
hydrogen atoms to fit the geometry.
If an atom only has one connection (terminal), it calculates bond distance and looks for
that bond distance in the bond_table and assigns a hybridization based on geometry.
Send any questions or comments to:   dmkilburg@gmail.com

Modified by John D. Chodera at MSKCC to utilize the OpenEye toolkit to enumerate protonation/tautomer
states and assign canonical AM1-BCC charges, and to compute pKas with epik from the Schrodinger Suite.
email: john.chodera@choderalab.org

Examples
--------

Extract bosutinib (resname DB8 from chain A) from Abl (PDB code 3UE4) and parameterize it.

> $SCHRODINGER/run tpl_maker_am1bcc.py -p rcsb-pdbs/3UE4.pdb DB8 A 7

Extract imatinib (resname STI from chain A) from Abl (PDB code 2HYY) and parameterize it.

> $SCHRODINGER/run tpl_maker_am1bcc.py -p rcsb-pdbs/2HYY.pdb STI A 7

"""

import os
import sys
import math
import datetime
import commands

from openeye import oechem, oequacpac, oeomega  # Requires OpenEye toolkit

class Atom(object):
    def __init__(self,name,idnum,element):
        self.name = name.strip()
        self.coords = []
        self.idnum = idnum
        self.connects = []
        self.element = element.strip()
        self.hybrid = ""

class Conformer(object):
    """
    Object representing a conformer with an associated pKa.

    """

    def __init__(self, name, charge, molecule, pKa=None, state_penalty=None):
        """
        Create a Conformer object.

        Parameters
        ----------
        name : str
            The short name of this conformer variant (usually three-letter ligand code + number)
        charge : int
            The corresponding integral total charge of this conformer.
        molecule : openeye.oechem.OEMol
            Molecule corresponding to this conformer.
        pKa : float, optional, default=None
            The corresponding pKa of this conformer.

        """
        self.name = name
        self.charge = charge
        self.molecule = molecule
        self.pKa = pKa
        self.state_penalty = state_penalty

        # Determine number of atoms
        atoms = [ atom for atom in molecule.GetAtoms() ]
        self.natoms = len(atoms)

class Pdb(object):
    """
    Storage for PDB-derived ligand information.

    """

    def add_atom(self,atom):
        if atom.idnum > self.max_idnum:
            self.max_idnum = atom.idnum
        self.atombyid[atom.idnum] = atom
        self.atom_list.append(atom)

    def add_hetatm(self,line,options):
        if options.ligand == line[17:20] and options.chain_identifier == line[21]:
            self.lines.append(line) # store this PDB record, since it belongs to the ligand
            atom = Atom(line[12:16],int(line[6:11].strip()),line[76:78])
            atom.coords.append(float(line[30:38]))
            atom.coords.append(float(line[38:46]))
            atom.coords.append(float(line[46:54]))
            self.add_atom(atom)

    def add_hydrogen(self,connect):
        self.max_idnum += 1
        atom = Atom("H"+str(self.hydrogen_count),self.max_idnum,'H')
        self.hydrogen_count += 1
        atom.connects.append(connect)
        atom.hybrid = 's'
        self.add_atom(atom)
        return atom

    def add_conect(self,line):
        idnum = int(line[6:11])
        if self.atombyid.has_key(idnum):
            self.lines.append(line) # store this PDB record, since it belongs to the ligand
            atom = self.atombyid[idnum]
            for i in [11,16,21,26]:
                c = line[i:i + 5].strip()
                if c != '':
                    connect_id = int(c)
                    if self.atombyid.has_key(connect_id):
                        atom.connects.append(connect_id)
        return

    def read_pdb(self,options,file):
        for line in file:
            if line.startswith("CONECT"):
                self.add_conect(line)
            elif line.startswith("HETATM"):
                self.add_hetatm(line,options)

    def read_ideal_pdb(self,options,file):
        for line in file:
            if line.startswith("ATOM  "):
                self.add_hetatm(line,options)
            elif line.startswith("CONECT"):
                self.add_conect(line)

    def read(self,options):
        with open(options.filename, 'r') as file:
            if options.ideal:
                self.read_ideal_pdb(options,file)
            else:
                self.read_pdb(options,file)

    def __init__(self,options):
        self.lines = list() # Collection of HETATM and CONECT records from PDB file for ligand.
        self.atom_list = []
        self.atombyid = {}
        self.max_idnum = 0
        self.hydrogen_count = 1
        self.read(options)
        self.pdb_extract = ''.join(self.lines) # concatenatd lines from PDB

        # DEBUG
        print "options.ligand = '%s'" % options.ligand
        print "options.chain_identifier = '%s'" % options.chain_identifier
        print self.pdb_extract

################ Function Definitions ##################

def write_comment_header(options,tpl):
    tpl.write('####################################\n')
    tpl.write('# Topology File for:\n')
    tpl.write('# {}\n'.format(options.ligand))
    tpl.write('# Extracted from: {}\n'.format(options.filename))
    tpl.write('#\n')
    tpl.write('# Created on: {}\n'.format(datetime.date.today()))
    tpl.write('#\n')
    tpl.write('# Created with: tpl_maker_am1bcc.py\n')
    tpl.write('####################################\n')
    tpl.write('\n')
    return

def assign_canonical_am1bcc_charges(molecule):
    """
    Assign canonical AM1-BCC charges to molecule.

    From canonical AM1-BCC recipe:
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        Molecule is modified in place.

    Raises
    ------
    An Exception is raised if no charges were set due to a failure.

    """

    # Create temporary copy of molecule to parameterize.
    expanded_molecule = oechem.OEMol(molecule)

    # Initialize Omega.
    omega = oeomega.OEOmega()
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)
    omega.SetSampleHydrogens(True)
    eWindow = 15.0
    omega.SetEnergyWindow(eWindow)
    omega.SetMaxConfs(800)
    omega.SetRMSThreshold(1.0)
    if not omega(expanded_molecule):
        raise ValueError("Cannot expand molecule with Omega Toolkit!")

    oequacpac.OEAssignPartialCharges(expanded_molecule, oequacpac.OECharges_AM1BCCSym)

    # Check whether charges were set.
    total_charge = 0.0
    nonzero_charges_detected = False
    for (src_atom, dest_atom) in zip(expanded_molecule.GetAtoms(), molecule.GetAtoms()):
        if src_atom.GetPartialCharge() != dest_atom.GetPartialCharge():
            nonzero_charges_detected = True
        total_charge += src_atom.GetPartialCharge()
    if not nonzero_charges_detected:
        raise ValueError("Charge assignment failure.")

    print " Total charge: %12.8f" % total_charge

    # Copy charges back to original molecule.
    for (src_atom, dest_atom) in zip(expanded_molecule.GetAtoms(), molecule.GetAtoms()):
        dest_atom.SetFormalCharge(src_atom.GetFormalCharge())
        dest_atom.SetPartialCharge(src_atom.GetPartialCharge())

    return

def assign_simple_am1bcc_charges(molecule, verbose=True):
    """
    Assign AM1-BCC charges to molecule using a single conformer, rather than the canonical scheme.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        Molecule is modified in place.

    Raises
    ------
    An Exception is raised if no charges were set due to a failure.

    """

    # Create temporary copy of molecule to parameterize.
    expanded_molecule = oechem.OEMol(molecule)

    oequacpac.OEAssignPartialCharges(expanded_molecule, oequacpac.OECharges_AM1BCC)

    # Check whether charges were set.
    nonzero_charges_detected = False
    total_charge = 0.0
    for (src_atom, dest_atom) in zip(expanded_molecule.GetAtoms(), molecule.GetAtoms()):
        if src_atom.GetPartialCharge() != dest_atom.GetPartialCharge():
            nonzero_charges_detected = True
        total_charge += src_atom.GetPartialCharge()
    if not nonzero_charges_detected:
        raise Exception("Charge assignment failure.")
    if verbose: print " Total charge: %12.8f" % total_charge

    # Copy charges back to original molecule.
    for (src_atom, dest_atom) in zip(expanded_molecule.GetAtoms(), molecule.GetAtoms()):
        dest_atom.SetPartialCharge(src_atom.GetPartialCharge())
        dest_atom.SetFormalCharge(src_atom.GetFormalCharge())

    return

def mk_single_conformer_epik(mol_name, molecule, verbose=True, pH=7):
    """
    Enumerate the list of conformers and associated properties for each protonation and
    tautomeric state using epik from the Schrodinger Suite, and save the conformer with
    the lowest penalty.

    Parameters
    ----------
    mol_name : string
        Molecule's name.
    molecule : openeye.oechem
        The molecule whose protomer and tautomer states are to be enumerated.
    pH : float, optional, default=7.0
        pH to use for conformer enumeration

    Returns
    -------
    molecule : openeye.oechem.OEMol
        The conformer with the lowest penalty.

    """

    # Write mol2 file.
    if verbose:
        print "Writing input file as mol2..."
    outmol = oechem.OEMol(molecule)
    ofs = oechem.oemolostream()
    ofs.open('epik-input.mol2')
    oechem.OEWriteMolecule(ofs, outmol)
    ofs.close()

    # Run Schrodinger's epik
    try:
        import run_epik
        run_epik.main(verbose, pH)
    except ImportError:
        import subprocess
        schrodinger_run = os.path.join(os.environ['SCHRODINGER'], 'run')
        run_epik_script = os.path.join(os.path.dirname(__file__), 'run_epik.py')
        cmd = schrodinger_run + ' ' + run_epik_script + ' {!r} {:f}'.format(verbose, pH)
        return_code = subprocess.call(cmd, shell=True)
        if return_code != 0:
            raise RuntimeError('{!r} failed, exit status: {:d}'.format(cmd, return_code))

    # Find minimum charge.
    ifs = oechem.oemolistream()
    ifs.open('epik-output.mol2')
    molecule = oechem.OEGraphMol()
    min_formal_charge = 1000
    while oechem.OEReadMolecule(ifs, molecule):
        # Check aromaticity.
        oechem.OEAssignAromaticFlags(molecule)

        # Assign formal charge
        oechem.OEAssignFormalCharges(molecule)
        formal_charge = 0
        for atom in molecule.GetAtoms():
            formal_charge += atom.GetFormalCharge()

        # Keep most negative formal charge
        min_formal_charge = min(min_formal_charge, formal_charge)

    ifs.close()
    if verbose:
        print "Minimum formal charge = %d" % min_formal_charge

    # Read conformer with lowest penalty from SDF (converted from Epik).
    if verbose:
        print "Reading conformers from SDF..."
    ifs_sdf = oechem.oemolistream()
    ifs_sdf.SetFormat(oechem.OEFormat_SDF)
    ifs_sdf.open('epik-output.sdf')
    sdf_molecule = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs_sdf, sdf_molecule)  # Read the molecule with lowest penalty
    ifs_sdf.close()

    # Get Epik data.
    epik_Ionization_Penalty = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty"))
    epik_Ionization_Penalty_Charging = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty_Charging"))
    epik_Ionization_Penalty_Neutral = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty_Neutral"))
    epik_State_Penalty = float(oechem.OEGetSDData(sdf_molecule, "r_epik_State_Penalty"))
    epik_Tot_Q = int(oechem.OEGetSDData(sdf_molecule, "i_epik_Tot_Q"))

    # Read conformer with lowest penalty from MOL2 (converted from Epik).
    ifs_mol2 = oechem.oemolistream()
    ifs_mol2.open('epik-output.mol2')
    mol2_molecule = oechem.OEMol()
    oechem.OEReadMolecule(ifs_mol2, mol2_molecule)
    ifs_mol2.close()

    # Make a copy of the mol2 molecule.
    molecule = oechem.OEMol(mol2_molecule)
    oechem.OEAssignAromaticFlags(molecule)
    molecule.SetTitle(mol_name)

    # Assign formal charge
    oechem.OEAssignFormalCharges(molecule)
    formal_charge = 0.0
    for atom in molecule.GetAtoms():
        formal_charge += atom.GetFormalCharge()
    if verbose:
        print "formal charge: %d" % formal_charge

    # Assign canonical AM1BCC charges.
    if verbose:
        print "Assigning AM1-BCC charges..."
    assign_canonical_am1bcc_charges(molecule)

    # Save result in mol2 file
    if verbose:
        print "Writing %s to mol2..." % mol_name
    ofs = oechem.oemolostream()
    ofs.open(mol_name + '.mol2')
    oechem.OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Compute number of protons.
    nprotons = epik_Tot_Q - min_formal_charge + 1

    # Compute effective pKa.
    import numpy as np
    kT = 298 * 6.022e23 * 1.381e-23 / 4184 # kcal/mol for 298 K
    pKa = pH - epik_State_Penalty / (nprotons * kT * np.log(10))
    print "effective pKa = %8.3f" % pKa

    # DEBUG
    print "%24s : pKa penalty %8.3f kcal/mol | tautomer penalty %8.3f kcal/mol | total state penalty %8.3f" % (mol_name, epik_Ionization_Penalty, epik_State_Penalty - epik_Ionization_Penalty, epik_State_Penalty)
    print epik_Tot_Q

    return molecule


def mk_conformers_epik(options, molecule, maxconf=99, verbose=True, pH=7):
    """
    Enumerate the list of conformers and associated properties for each protonation and tautomeric state using epik from the Schrodinger Suite.
    Parameters
    ----------
    options
    molecule : openeye.oechem
        The molecule read from the PDB whose protomer and tautomer states are to be enumerated.
    maxconf : int, optional, default=128
        Maximum number of protomers/tautomers to generate.
    pH : float, optional, default=7.0
        pH to use for conformer enumeration
    Returns
    -------
    conformers : list of Conformer
        The list of protomers/tautomers generated.
    """

    from schrodinger import structure # Requires Schrodinger Suite

    # Write mol2 file.
    if verbose: print "Writing input file as mol2..."
    outmol = oechem.OEMol(molecule)
    ofs = oechem.oemolostream()
    ofs.open('epik-input.mol2')
    oechem.OEWriteMolecule(ofs, outmol)
    ofs.close()
    # Use low level writer to get atom names correct.
    ofs = oechem.oemolostream()
    ofs.open('epik-input.mol2')
    for (dest_atom, src_atom) in zip(outmol.GetAtoms(), molecule.GetAtoms()):
        dest_atom.SetName(src_atom.GetName())
    oechem.OEWriteMol2File(ofs, outmol, True)
    ofs.close()

    # Write mol2 file.
    if verbose: print "Writing input file as sdf..."
    outmol = oechem.OEMol(molecule)
    ofs = oechem.oemolostream()
    ofs.open('epik-input.sdf')
    oechem.OEWriteMolecule(ofs, outmol)
    ofs.close()

    # Write pdb file.
    if verbose: print "Writing input file as pdb..."
    outmol = oechem.OEMol(molecule)
    ofs = oechem.oemolostream()
    ofs.open('epik-input.pdb')
    oechem.OEWriteMolecule(ofs, outmol)
    ofs.close()

    # Write input for epik.
    if verbose: print "Converting input file to Maestro format..."
    reader = structure.StructureReader("epik-input.mol2")
    writer = structure.StructureWriter("epik-input.mae")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()

    # Run epik to enumerate protomers/tautomers and get associated state penalties.
    if verbose: print "Running Epik..."
    cmd = '%s/epik -imae epik-input.mae -omae epik-output.mae -pht 10.0 -ms 100 -nt -pKa_atom -ph %f -WAIT' % (os.environ['SCHRODINGER'], pH)
    output = commands.getoutput(cmd)
    if verbose: print output

    # Convert output from epik from .mae to .sdf.
    if verbose: print "Converting output file to SDF..."
    reader = structure.StructureReader("epik-output.mae")
    writer = structure.StructureWriter("epik-output.sdf")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()

    # Also convert to .mol2.
    if verbose: print "Converting output file to MOL2..."
    reader = structure.StructureReader("epik-output.mae")
    writer = structure.StructureWriter("epik-output.mol2")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()

    # Find minimum charge.
    ifs = oechem.oemolistream()
    ifs.open('epik-output.mol2')
    molecule = oechem.OEGraphMol()
    min_formal_charge = 1000
    while oechem.OEReadMolecule(ifs, molecule):
        # Check aromaticity.
        oechem.OEAssignAromaticFlags(molecule)

        # Assign formal charge
        oechem.OEAssignFormalCharges(molecule)
        formal_charge = 0
        for atom in molecule.GetAtoms():
            formal_charge += atom.GetFormalCharge()
        # Keep most negative formal charge
        min_formal_charge = min(min_formal_charge, formal_charge)

    ifs.close()
    if verbose: print "Minimum formal charge = %d" % min_formal_charge

    # Read conformers from SDF and mol2 (converted from Epik).
    if verbose: print "Reading conformers from SDF..."
    ifs_sdf = oechem.oemolistream()
    ifs_sdf.SetFormat(oechem.OEFormat_SDF)
    ifs_sdf.open('epik-output.sdf')
    sdf_molecule = oechem.OEGraphMol()

    ifs_mol2 = oechem.oemolistream()
    ifs_mol2.open('epik-output.mol2')
    mol2_molecule = oechem.OEGraphMol()

    conformer_index = 1
    conformers = list()
    while oechem.OEReadMolecule(ifs_sdf, sdf_molecule):
        if verbose: print "Conformer %d" % conformer_index

        # Read corresponding mol2 molecule.
        oechem.OEReadMolecule(ifs_mol2, mol2_molecule)
        oechem.OEAssignAromaticFlags(mol2_molecule) # check aromaticity

        # Make a copy of the mol2 molecule.
        molecule = oechem.OEMol(mol2_molecule)

        # Set name
        name = options.ligand+'%02d' % conformer_index
        molecule.SetTitle(name)

        # Assign formal charge
        oechem.OEAssignFormalCharges(molecule)
        formal_charge = 0.0
        for atom in molecule.GetAtoms():
            formal_charge += atom.GetFormalCharge()
        if verbose: print "formal charge: %d" % formal_charge

        # DEBUG: Write mol2 file before assigning charges.
        if verbose: print "Writing %s to mol2..." % name
        outmol = oechem.OEMol(molecule)
        ofs = oechem.oemolostream()
        ofs.open(name + '.mol2')
        oechem.OEWriteMolecule(ofs, outmol)
        ofs.close()

        # Assign canonical AM1BCC charges.
        try:
            if verbose: print "Assigning AM1-BCC charges..."
            #assign_canonical_am1bcc_charges(molecule)
            assign_simple_am1bcc_charges(molecule)
        except Exception as e:
            print str(e)
            continue

        # Get Epik data.
        epik_Ionization_Penalty = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty"))
        epik_Ionization_Penalty_Charging = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty_Charging"))
        epik_Ionization_Penalty_Neutral = float(oechem.OEGetSDData(sdf_molecule, "r_epik_Ionization_Penalty_Neutral"))
        epik_State_Penalty = float(oechem.OEGetSDData(sdf_molecule, "r_epik_State_Penalty"))
        epik_Tot_Q = int(oechem.OEGetSDData(sdf_molecule, "i_epik_Tot_Q"))

        # Compute number of protons.
        nprotons = epik_Tot_Q - min_formal_charge + 1

        # Compute effective pKa.
        import numpy as np
        kT = 298 * 6.022e23 * 1.381e-23 / 4184 # kcal/mol for 298 K
        pKa = options.pH - epik_State_Penalty / (nprotons * kT * np.log(10))
        print "effective pKa = %8.3f" % pKa

        # DEBUG
        print "%24s : pKa penalty %8.3f kcal/mol | tautomer penalty %8.3f kcal/mol | total state penalty %8.3f\n" % (name, epik_Ionization_Penalty, epik_State_Penalty - epik_Ionization_Penalty, epik_State_Penalty)

        # Create a conformer and append it to the list.
        conformer = Conformer(name, epik_Tot_Q, molecule, state_penalty=epik_State_Penalty)
        conformers.append(conformer)
        print epik_Tot_Q # DEBUG
        # Increment counter.
        conformer_index += 1

    ifs_sdf.close()
    ifs_mol2.close()

    if verbose: print "%d protomer/tautomer states were enumerated" % len(conformers)

    return conformers


def mk_conformers(options, molecule, maxconf=99, verbose=True):
    """
    Enumerate the list of conformers and associated properties for each protonation and tautomeric state.

    Parameters
    ----------
    options
    molecule : openeye.oechem
        The molecule read from the PDB whose protomer and tautomer states are to be enumerated.
    maxconf : int, optional, default=128
        Maximum number of protomers/tautomers to generate.

    Returns
    -------
    conformers : list of Conformer
        The list of protomers/tautomers generated.

    """
    conformers = list()

    formalChargeOptions = oequacpac.OEFormalChargeOptions(maxconf)
    tautomerOptions = oequacpac.OETautomerOptions(maxconf)
    # enumerate pka states
    index = 0
    for pkaState in oequacpac.OEEnumerateFormalCharges(molecule, formalChargeOptions):
        if (index > maxconf): break
        # enumerate tautomers of pka states
        for tautomer in oequacpac.OEEnumerateTautomers(pkaState, tautomerOptions):
            if (index > maxconf): break
            # Assign conformer name.
            name = options.ligand+'%02d' % index
            tautomer.SetTitle(name)
            # DEBUG: Write molecule before charging.
            ofs = oechem.oemolostream()
            ofs.open(name + '-before.mol2')
            outmol = oechem.OEMol(tautomer)
            oechem.OEWriteMolecule(ofs, outmol)
            ofs.close()
            # Compute formal charge.
            oechem.OEAssignFormalCharges(tautomer)
            formal_charge = 0.0
            for atom in tautomer.GetAtoms():
                formal_charge += atom.GetFormalCharge()
            if verbose: print "formal charge: %d" % formal_charge
            # Assign canonical AM1BCC charges.
            assign_canonical_am1bcc_charges(tautomer)
            # DEBUG: Write molecule after charging.
            ofs = oechem.oemolostream()
            outmol = oechem.OEMol(tautomer)
            ofs.open(name + '-after.mol2')
            oechem.OEWriteMolecule(ofs, outmol)
            ofs.close()
            # Create conformer.
            conformer = Conformer(name, formal_charge, tautomer)
            # Append to protomer/tautomer list.
            conformers.append(conformer)
            # Keep count of protomers/tautomers.
            index += 1
            # Show output if requested.
            if verbose:
                print "%12s %+5d %8d atoms" % (name, conformer.charge, conformer.natoms)

    if verbose: print "%d protomer/tautomer states were enumerated" % len(conformers)

    return conformers

def write_conformers(options,tpl,conformers):
    tpl.write('CONFLIST {}        {}BK '.format(options.ligand,options.ligand))
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write('{} '.format(conformer.name))
    tpl.write('{} '.format(options.ligand + "DM"))
    tpl.write('\n')
    tpl.write('\n')
    return

def mk_vectors(atom1,atom2,atom3):
    vect1 = []
    vect2 = []
    for i in range(3):
        vect1.append(float(atom2.coords[i])-float(atom1.coords[i]))
        vect2.append(float(atom3.coords[i])-float(atom1.coords[i]))
    return vect1,vect2

def mk_dot_prod(vect1,vect2):
    dot_prod = 0
    for i in range(3):
        dot_prod = dot_prod + vect1[i]*vect2[i]
    return dot_prod

def magnitude(vector):
    mag = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    return mag

def find_angle(vect1,vect2):
    try:
        cos_ang = mk_dot_prod(vect1,vect2)/(magnitude(vect1)*magnitude(vect2))
        angle = math.degrees(math.acos(cos_ang))
        return angle
    except:
        return None

## Bond distances from Allen, F.H.; Kennard, O. et al. (1987) Tables of 
## Bond Lengths determined by X-Ray and Neutron Diffraction. Part I . 
## Bond Lengths in Organic Compounds. J. CHEM. SOC. PERKIN TRANS

## [Atom 1, Atom 2, lower limit(bond length), upper limit(bond length), hybridization of atom 1]

bond_table = [
    [ "C", "C", 1.514, 1.588, 'sp3'], ## C-C
    [ "C", "C", 1.501, 1.513, 'sp2'], ## C=C
    [ "C", "C", 1.465, 1.473, 'sp'],  ## C-tb-C, tb = triple bond
    [ "C", "N", 1.487, 1.510, 'sp3'], ## C-N
    [ "C", "N", 1.278, 1.432, 'sp2'], ## C=N
    [ "C", "N", 1.134, 1.146, 'sp' ], ## C-tb-N
    [ "C", "O", 1.389, 1.50, 'sp3'],  ## C-O
    [ "C", "O", 1.180, 1.387, 'sp2'], ## C=O
    [ "C", "S", 1.744, 1.865, 'sp3'], ## C-S
    [ "C", "S", 1.710, 1.743, 'sp2'], ## C=S
    [ "C", "S", 1.629, 1.679, 'sp'],  ## C-tb-S
    [ "C", "P", 0.00, 2.00, 'sp3'],
    [ "C", "BR", 0.00, 3.00,'sp3'],
    [ "C", "CL", 0.00, 3.00, 'sp3'],
    [ "C", "F", 0.00 , 3.00, 'sp3'],
    [ "C", "I", 0.00, 3.00, 'sp3'],
    [ "O", "C", 1.30, 3.00, 'sp3'],
    [ "O", "C", 0.00, 1.29, 'sp2'],
    [ "O", "S", 1.476, 3.00, 'sp3'],
    [ "O", "S", 0.00, 1.476, 'sp2'],
    [ "O", "P", 1.508, 3.00, 'sp3'],
    [ "O", "P", 0.00, 1.508, 'sp2'],
    [ "O", "N", 1.30, 3.00, 'sp3'],
    [ "O", "N", 0.00, 1.29, 'sp2'],
    [ "O", "O", 0.00, 3.00, 'sp3'],
    [ "N", "C", 1.468, 1.511, 'sp3'],
    [ "N", "C", 1.27, 1.432, 'sp2'],
    [ "N", "C", 1.13, 1.16, 'sp'],
    [ "N", "N", 1.300, 1.455, 'sp3'],
    [ "N", "N", 1.122, 1.299, 'sp2'],
    [ "N", "O", 1.24, 1.47, 'sp3'],
    [ "N", "O", 1.15, 1.239, 'sp2'],
    [ "N", "P", 1.625, 1.70, 'sp3'],
    [ "N", "P", 1.560, 1.61, 'sp2'],
    [ "N", "S", 1.59, 1.75, 'sp3'],
    [ "N", "S", 1.52, 1.58, 'sp2'],
    [ "N", "SI", 0.00, 3.00, 'sp3'],
    [ "P", "P", 2.00, 3.00, '???'],
    [ "P", "P", 0.00, 1.99, 'sp3'],
    [ "P", "N", 1.560, 1.61, 'sp3'],
    [ "P", "C", 0.00, 3.00, 'sp3'],
    [ "P", "O", 0.00, 3.00, 'sp3'],
    [ "P", "S", 0.00, 3.00, 'sp3'],
    [ "S", "N", 1.50, 1.6, 'sp'],
    [ "S", "N", 1.6, 3.00, 'sp2'],
    [ "S", "O", 1.50, 1.60, 'sp2'],
    [ "S", "O", 1.65, 3.00, 'sp'],
    [ "S", "N", 1.50, 3.00, 'sp2'],
    [ "S", "N", 0.00, 1.45, 'sp'],
    [ "S", "C", 0.00, 3.00, 'sp2']
]
def search_bond_table(a1,a2,distance):
    for row in bond_table:
        if row[0] == a1.element and row[1] == a2.element \
                and row[2] <= distance and distance <= row[3]:
            return row[4]
    return 'unknownB'

## for use when you do not have to add hydrogens

hydrogen_table2 = [
    ['C', 4, 'sp3'],
    ['C', 3, 'sp2'],
    ['C', 2, 'sp'],
    ['N', 3, 'sp3'],
    ['N', 2, 'sp2'],
    ['O', 2, 'sp3'],
    ['O', 1, 'sp2'],
    ['H', 1, 's'],
    ['S', 4, 'sp3'],
    ['P', 4, 'sp3']
]

def search_hyd2_table(a1,connects):
    for row in hydrogen_table2:
        if row[0] == a1 and row[1] == connects:
            return row[2]
    return 'unknownH'

def add_hybrids(pdb):
    for atom in pdb.atom_list:
        connects = len(atom.connects)
        atom.hybrid = search_hyd2_table(atom.element,connects)
    return
    
## [Atom, number of atoms connected, hybridization, number of hydrogens to add]
hydrogen_table = [
    ['C', 4, 'sp3', 0],
    ['C', 3, 'sp3', 1],
    ['C', 3, 'sp2', 0],
    ['C', 2, 'sp3', 2],
    ['C', 2, 'sp2', 1],
    ['C', 2, 'sp', 0],
    ['C', 1, 'sp3', 3],
    ['C', 1, 'sp2', 2],
    ['C', 1, 'sp', 1],
    ['N', 3, 'sp3', 0],
    ['N', 3, 'sp2', 0],
    ['N', 2, 'sp3', 1],
    ['N', 2, 'sp2', 1],
    ['N', 2, 'sp', 0],
    ['N', 1, 'sp3', 2],
    ['N', 1, 'sp2', 1],
    ['N', 1, 'sp', 0],
    ['O', 2, 'sp3', 0],
    ['O', 1, 'sp3', 1],
    ['O', 1, 'sp2', 0],
    ['S', 4, 'sp3', 0],
    ['S', 3, 'sp3', 1],
    ['S', 3, 'sp2', 0],
    ['S', 2, 'sp3', 2],
    ['S', 2, 'sp2', 1],
    ['S', 1, 'sp2', 1],
    ['P', 4, 'sp3', 0],
    ['P', 3, 'sp3', 1],
    ['P', 3, 'sp2', 0],
    ['P', 2, 'sp3', 2],
    ['P', 2, 'sp2', 1],
    ['P', 2, 'sp', 0],
    ['P', 1, 'sp3', 3],
    ['B', 3, 'sp3', 0],
    ['B', 3, 'sp2', 0],
    ['B', 2, 'sp3', 1],
    ['B', 2, 'sp2', 1],
]

def search_hydrogen_table(atom):
    for row in hydrogen_table:
        if row[0] == atom.element and row[2] == atom.hybrid \
                and row[1] == len(atom.connects):
            return row[3]
    return 0

def add_hydrogens(pdb):
    for atom1 in pdb.atom_list:
        con_len = len(atom1.connects)
        if con_len >= 2:
            atom2 = pdb.atombyid[atom1.connects[0]]
            atom3 = pdb.atombyid[atom1.connects[1]]
            (vector1, vector2) = mk_vectors(atom1,atom2,atom3)
            angle = find_angle(vector1,vector2)
            atom1.hybrid = geometry(angle)
        elif con_len == 1 and atom1.element != 'H':
            atom2 = pdb.atombyid[atom1.connects[0]]
            bond_len = bond_length(atom1,atom2)
            orbital = search_bond_table(atom1,atom2,bond_len)
            atom1.hybrid = orbital

        num_hydrogens = search_hydrogen_table(atom1)
        for i in range(num_hydrogens):
            h_atom = pdb.add_hydrogen(atom1.idnum)
            atom1.connects.append(h_atom.idnum)

def bond_length(atom1,atom2):
    vector = []
    for i in range(3):
        vector.append(atom1.coords[i] - atom2.coords[i])
    length = magnitude(vector)
    return length

# Returns hybridization based on angles.
def geometry(angle):
    if angle > 106 and angle < 115:
        return 'sp3'
    if angle < 183 and angle > 177:
        return 'sp'
    if angle < 125 and angle > 115:
        return 'sp2'
    else:
        return 'unknownA'


def write_natom(options,tpl,pdb,conformers):
    """
    Write block describing number of atoms present in each conformer.

    Parameters
    ----------
    options :
        Options object.
    tpl : file
        tpl file handle to write to.
    pdb : Pdb
        PDB object for input ligand.
    conformers : list of Conformer
        Conformers to write information for.

    """
    template = '{0:9}{1:11}{2:1}\n'
    tpl.write(template.format('NATOM',options.ligand+"BK", '0'))

    if options.reverse_order:
        for conformer in reversed(conformers):
            tpl.write(template.format('NATOM',
                                      conformer.name, conformer.natoms))
    else:
        for conformer in conformers:
            tpl.write(template.format('NATOM',
                                      conformer.name, conformer.natoms))

    tpl.write(template.format('NATOM',options.ligand+"DM", '0'))
    tpl.write('\n')

    return

def o_connected(atom,pdb):
    for idnum in atom.connects:
        if pdb.atombyid[idnum].element == 'O':
            return True
    return False

def write_atoms(options,tpl,pdb,conformers,printer):
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        count = 0
        for atom in conformer.molecule.GetAtoms():
            printer(tpl,conformer,atom,count)
            count += 1
        tpl.write('\n')
    return

def write_iatom(options,tpl,pdb,conformers):
    def printer(tpl,conformer,atom,count):
        template = '{0:9s}{1:6s}{2:>5s}{3:5d}\n'
        tpl.write(template.format('IATOM', conformer.name,
                                  atom.GetName(), count))
        return
    write_atoms(options,tpl,pdb,conformers,printer)
    return

def write_atomname(options,tpl,pdb,conformers):
    def printer(tpl,conformer,atom,count):
        #template = '{0:9}{1:8}{2:>2}  {3:>5}\n'
        template = '{0:9s}{1:6s}{2:5d}{3:>4s}\n'
        tpl.write(template.format('ATOMNAME', conformer.name, \
                                  count,atom.GetName()))
        return
    write_atoms(options,tpl,pdb,conformers,printer)
    return

def write_atom_param_section(options,tpl,pdb,conformers,vdw_dict):
    tpl.write("# Atom Parameters:\n")
    tpl.write("# Van Der Waals Radii. See source for reference\n")
    def printer(tpl,conformer,atom,count):
        #template = '{0:9}{1:7}{2:5} {3:7}\n'
        template = '{0:9s}{1:6s}{2:>5s}{3:+012.8f}\n'
        element = oechem.OEGetAtomicSymbol(atom.GetAtomicNum()).upper()
        tpl.write(template.format("RADIUS",conformer.name, \
                                  atom.GetName(),vdw_dict[element]))
        return
    write_atoms(options,tpl,pdb,conformers,printer)
    return

def write_charges(options,tpl,pdb,conformers):
    def printer(tpl,conformer,atom,count):
        template = '{0:9s}{1:6s}{2:>5s}{3:+012.8f}\n'
        charge = atom.GetPartialCharge()
        tpl.write(template.format("CHARGE",conformer.name, \
                                  atom.GetName(), charge))
        return
    write_atoms(options,tpl,pdb,conformers,printer)
    return

def write_sect1_header(tpl):
    tpl.write('# 1. Basic Conformer Information:\n')
    tpl.write('# Number of protons and electrons, pKa, Em, and Reaction Field Energy (RXN)\n')
    tpl.write('')
    return

def write_proton(options,tpl,conformers):
    tpl.write('# PROTON SECTION: PROTON means charge:\n')
    template = '{0:9}{1:11}{2:5}\n'
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write(template.format("PROTON",conformer.name, \
                                      str(conformer.charge)))
    tpl.write(template.format("PROTON",options.ligand+"DM", '0'))
    tpl.write('\n')
    return

def write_pka(options,tpl,conformers):
    tpl.write('# Solution pKa Section: pKa data from CRC Handbook of Chemistry and Physics\n')
    template = '{0:9s}{1:6s}     {2:8.3f}\n'
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        if conformer.pKa == None:
            tpl.write(template.format("PKA", conformer.name, options.pH))
        else:
            tpl.write(template.format("PKA", conformer.name, conformer.pKa))
    tpl.write(template.format("PKA", options.ligand+"DM", 0.0))
    tpl.write('\n')
    return

def write_electron(options,tpl,conformers):
    tpl.write("#ELECTRON SECTION:\n")
    template = '{0:9}{1:11}{2:5}\n'
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write(template.format("ELECTRON",conformer.name,"0.0"))
    tpl.write(template.format("ELECTRON",options.ligand+"DM","0.0"))
    tpl.write('\n')
    return

def write_EM(options,tpl,conformers): 
    template = '{0:9}{1:11}{2:5}\n'
    tpl.write("# EM SECTION:\n")
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write(template.format("EM",conformer.name,"0.0"))
    tpl.write(template.format("EM",options.ligand+"DM","0.0"))
    tpl.write('\n')
    return

def write_RXN(tpl): 
    tpl.write("# REACTION FIELD ENERGY SECTION:\n")
    tpl.write('\n')
    return

def write_con_header(tpl): 
    string = "ires" + " " + "conn "
    tpl.write("#ONNECT" + "   " + "conf" + " " + "atom" + "  " + "orbital" \
                  + "  " +string*4 + "\n")
    tpl.write("#ONNECT" +" " + "|-----|----|---------|" + "----|"*10 + "\n")
    return

def write_con_section(options,tpl,pdb,conformers):
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write("#  " + conformer.name + "\n")
        write_con_header(tpl)
        template1 = '{0:9}{1:5} {2:^4} {3:^10}'
        template2 = '{0:^4} {1:^4} '
        for atom in conformer.molecule.GetAtoms():
            hyb_table = ['s', 'sp', 'sp2', 'sp3', 'sp3d', 'sp3d2']
            hybrid = hyb_table[oechem.OEGetHybridization(atom)]
            tpl.write(template1.format("CONNECT", conformer.name, atom.GetName(), hybrid))
            for bond in conformer.molecule.GetBonds():
                if bond.GetBgn() == atom:
                    tpl.write(template2.format("0", bond.GetEnd().GetName()))
                if bond.GetEnd() == atom:
                    tpl.write(template2.format("0", bond.GetBgn().GetName()))
            tpl.write('\n')
        tpl.write('\n')
    tpl.write('\n')
    return

import re

#################  Main  ###########################################

## vdw_dict contains vdw radius of atoms as written in
## Bondi,A. van Der Waals Volumes and Radii,
## J.Phys.Chem, 63, 3, 1964 and http://periodic.lanl.gov/index.shtml
## (Los Alamos National Lab Periodic Table)

vdw_dict = {'H':1.20,
            'He':1.40,
            'C':1.70,
            'N':1.55,
            'O':1.52,
            'F':1.47,
            'NE':1.54,
            'SI':2.10,
            'P':1.80,
            'S':1.80,
            'CL':1.75,
            'AR':1.88,
            'AS':1.85,
            'SE':1.90,
            'BR':1.85,
            'KR':2.02,
            'TE':2.06,
            'I':1.98,
            'XE':2.16,
            'ZN':1.39,
            'CU':1.40,
            'HG':1.55,
            'CD':1.58,
            'NI':1.63,
            'PD':1.63,
            'AU':1.66,
            'AG':1.72,
            'MG':1.73,
            'PT':1.75,
            'LI':1.82,
            'U':1.86,
            'GA':1.87,
            'PB':2.02,
            'SN':2.17,
            'NA':2.27,
            'FE':0.00,
            'K':2.75}

def create_openeye_molecule(pdb, options, verbose=True):
    """
    Create OpenEye molecule from PDB representation.

    The molecule will have hydrogens added and be normalized, but the overall geometry will not be altered.

    Parameters
    ----------
    pdb : Pdb
       The PDB-extracted entries for the ligand.

    Returns
    -------
    molecule : openeye.oechem.OEMol
        Molecule representation.
    options : options struct
        Options structure.

    """

    # Create a molecule container.
    molecule = oechem.OEGraphMol()

    # Open a PDB file reader from the stored PDB string representation of HETATM and CONECT records.
    print pdb.pdb_extract
    ifs = oechem.oemolistream()
    ifs.openstring(pdb.pdb_extract)
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_ALL
    ifs.SetFlavor(oechem.OEFormat_PDB, flavor)
    oechem.OEReadPDBFile(ifs, molecule)

    # Add explicit hydrogens.
    oechem.OEDetermineConnectivity(molecule)
    oechem.OEFindRingAtomsAndBonds(molecule)
    oechem.OEAssignAromaticFlags(molecule) # check aromaticity
    oechem.OEPerceiveBondOrders(molecule)

    # We must assign implicit hydrogens first so that the valence model will be correct.
    oechem.OEAssignImplicitHydrogens(molecule)
    oechem.OEAssignFormalCharges(molecule)

    # Now add explicit hydrogens.
    polarOnly = False
    set3D = True
    oechem.OEAddExplicitHydrogens(molecule, polarOnly, set3D)

    # TODO: Sequentially number hydrogen atoms.

    # Perceive stereochemostry.
    oechem.OEPerceiveChiral(molecule)

    # Set title.
    molecule.SetTitle(options.ligand)

    # Write out PDB form of this molecule.
    # TODO: Fix atom numbering.
    #if verbose: print "Writing input molecule as PDB..."
    #outmol = oechem.OEMol(molecule)
    #ofs = oechem.oemolostream()
    #flavor = oechem.OEOFlavor_Generic_Default | oechem.OEOFlavor_PDB_Default
    #ofs.SetFlavor(oechem.OEFormat_PDB, flavor)
    #ofs.open(options.ligand + '.pdb')
    #oechem.OEWriteMolecule(ofs, outmol)
    #ofs.close()

    # Write mol2 file for this molecule.
    if verbose: print "Writing input molecule as mol2..."
    outmol = oechem.OEMol(molecule)
    ofs = oechem.oemolostream()
    filename = options.ligand + '.mol2'
    ofs.open(filename)
    oechem.OEWriteMolecule(ofs, outmol)
    ofs.close()
    # Use low level writer to get atom names correct.
    ofs = oechem.oemolostream()
    ofs.open(filename)
    for (dest_atom, src_atom) in zip(outmol.GetAtoms(), molecule.GetAtoms()):
        dest_atom.SetName(src_atom.GetName())
    oechem.OEWriteMol2File(ofs, outmol, True)
    ofs.close()
    # Read and write in PDB format.
    if verbose: print "Converting mol2 to PDB..."
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    if ifs.open(options.ligand + '.mol2'):
        if ofs.open(options.ligand + '.pdb'):
            for mol in ifs.GetOEGraphMols():
                oechem.OEWriteMolecule(ofs, mol)

    return molecule

def write_extra(options,tpl,conformers):
    """
    """
#0123456789012345678901234567890123456789
#EXTRA    TYR-1      -0.80

    template = '{0:9s}{1:6s}{2:5s}{3:8.3f}\n'
    tpl.write("# EXTRA energy for tautomers:\n")
    if options.reverse_order:
        conformers = reversed(conformers)
    for conformer in conformers:
        tpl.write(template.format("EXTRA", conformer.name, "", conformer.state_penalty))
    tpl.write('\n')

    return

def write_tpl(options,tpl,pdb):
    # Create an OpenEye molecule from the PDB representation.
    molecule = create_openeye_molecule(pdb, options)

    # Write the MCCE2 .tpl file header.
    write_comment_header(options, tpl)

    # TODO: Figure out what this does.
    if options.ideal:
        add_hybrids(pdb)
    else:
        add_hydrogens(pdb)

    # Generate list of conformers with different protonation and tautomer states.
    #conformers = mk_conformers(options, molecule) # use OEProton protomers/tautomers
    conformers = mk_conformers_epik(options, molecule, pH=options.pH) # use Epik protomers/tautomers

    # Write the conformer definitions to the MCCE2 .tpl file.
    write_conformers(options,tpl,conformers)

    # Write the number of atoms associated with each conformer.
    write_natom(options,tpl,pdb,conformers)
    
    write_iatom(options,tpl,pdb,conformers)
    
    write_atomname(options,tpl,pdb,conformers)
    
    write_sect1_header(tpl)
    
    write_proton(options,tpl,conformers)
    
    write_pka(options,tpl,conformers)
    
    write_electron(options,tpl,conformers)
    
    write_EM(options,tpl,conformers)
    
    write_RXN(tpl)
    
    write_con_section(options,tpl,pdb,conformers)
    
    write_atom_param_section(options,tpl,pdb,conformers,vdw_dict)
    
    write_charges(options,tpl,pdb,conformers)

    write_extra(options,tpl,conformers)

    return

# Main

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', action="store", help = 'file to parse')
    parser.add_argument('ligand', action="store",
                        help = '3 letter code of ligand to extract from file')
    parser.add_argument('chain_identifier', action='store',
                        help = 'chain id for ligand')
    parser.add_argument('pH', action='store',type=float,
                        help = 'pH for which parameters are to be created.')
    parser.add_argument('-p', action='store_false',dest='ideal',default =True,
                        help = 'Extract ligand from PDB HETATM entries.')
    parser.add_argument('-r', action='store_false',dest='reverse_order', default=False,
                        help = 'print conformers starting from negative ions')
    options = parser.parse_args()

    with open(options.ligand+'.tpl','w') as tpl:
        write_tpl(options,tpl,Pdb(options))

if __name__ == "__main__":
    main ()
