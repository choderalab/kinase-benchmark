# Kinase benchmark

### Generate ligands structures
Run

```
python generate_ligands.py
```

The script takes all the SMILES in `clinical-kinase-inhibitors.csv` and uses the OpenEye Toolkit and Schrodinger’s epik to generate a structure and a protonation state for the ligands. The resulting mol2 files are saved in the folder `ligands`.

The environment variable `SCHRODINGER` must point to its installation directory (on Mac OS X it is usually `/opt/schrodinger/suites201X`).

This currently gives structures for all the FDA inhibitors except for Palbociclib which causes the following warnings

```
in function omega()
Warning: OE3DToAtomStereo is unable to perceive atom stereo from a flat geometry on atom 11 of molecule 'Palbociclib'
Warning: Palbociclib failed due to unspecified stereochemistry with strict stereo enabled

in function oequacpac.OEAssignPartialCharges()
Warning: BCC charge corrections will not be done for molecule Palbociclib
         because bcc parameters for C2--C3 (C2--C3) bond can't be assigned
```


### Setup yank simulation
You can setup Yank’s simulation for different kinase structures and inhibitors for example

```
python setup_yank.py --receptor abl/2HYY-pdbfixer.pdb --ligand Imatinib.mol2 --solvate
```

and

```
python setup_yank.py --receptor abl/3UE4-pdbfixer.pdb --ligand Bosutinib.mol2
```

create a setup folder to run Yank respectively on 2HYY and imatinib in explicit solvent and 3UE4 and bosutinib in implicit solvent.

The script uses a simple heuristic algorithm to ensure that the receptor and the ligand are very close and do not overlap at the beginning of the simulation.
