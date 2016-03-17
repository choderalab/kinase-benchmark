## Questions
* can repex overcome minor differences? (same conformation, different structures)
* what about different inactive conformations?
* what if the active/inactive kinetic barrier is too big?
* what about effect of modeled loops? (delete loops in resolved structure and remodel them)
* effect of SH2/SH3

## Combinations
* kinase targets: 2 * 3/4 * 1/2 * 2 * 2 = 24/64 combinations
  * Src, Abl: 2
  * which conformation? 3/4
    * source x-ray structs: sensitivity to initial structure, active and inactive
    * **(LATER)** MSM states: dissect contributions of reorganizing energy
  * which construct? 1/2
  * which protonation state? 2
    *DFG Asp (should not affect thermodynamics)
  * which phosphorylation state? 2
    *phosphotyrosine in activation loop (pTyr)
  * **(LATER)** mutants
* FDA approved inhibitors:
  * select noncovalent inhibitors in their active forms
  * which tautomeric/protonation state?
    *pick states with epik penalty <= 6kT
    *relevant epik conformation for imatinib is 0
  * **(LATER)** torsions parametrization?
* solvents:
  * explicit: PME, RF
  * implicit: GBSA-OBC (4-5 models)
  * ions: neutralize, set ionic strength
* force fields:
  * amber99sbildn, amber14
  * **(LATER)** charmm36
