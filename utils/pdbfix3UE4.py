from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

fixer = PDBFixer(pdbid='3UE4')

fixer.removeChains(chainIds=['B'])

# Without fixer.missingResidues = {}, fixer.addMissingAtoms() throw an exception
# and if I call fixer.findMissingResidues() several terminal residues are added
fixer.missingResidues = {}
fixer.findMissingAtoms()
fixer.addMissingAtoms()

fixer.removeHeterogens(keepWater=False)

#fixer.addMissingHydrogens(7.0)

PDBFile.writeFile(fixer.topology, fixer.positions, open('../kinases/abl/3UE4-pdbfixer.pdb', 'w'))
