# Utility Python functions

* pdbfix3UE4.py : used to extract the Abl domain saved in kinases/abl/.
* tpl_maker_am1bcc.py : this is based on [this script](https://github.com/choderalab/mcce-charges/blob/master/kinases/tpl_maker_am1bcc.py). I didn’t want to change anything so I just added a function (`mk_single_conformer_epik`) to create a single protonation state and assign canonical AM1-BCC charges.
* run_epik.py : isolate the part of the code that needs the shrodinger Python module so that the script can be run both with the normal python interpreter or Schrodinger’s.