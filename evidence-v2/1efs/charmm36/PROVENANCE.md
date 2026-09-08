# 1efs / CHARMM36

Structure: `examples/1efs_heteroduplex/1efs.pdb` at `ea00b023b9d853571a6cc9828e2cf1bb7098dc90` (blob `9716d2b2ee3609ae616dcd5cbd28a764cc8bf61e`), rebuilt as a 5'-OH/3'-OH strand for the CHARMM36 route: 816 atoms, 6 phosphate atoms removed, 4 hydroxyl H added (`topology/build_1efs_charmm36.log`). Deck `engine/params/db_charmm36_na` at `ea00b023b9d853571a6cc9828e2cf1bb7098dc90`. MOSAICS binary sha256 `8b6da340db0f7438bc9b1d195dcfba58db83b2140fc29bceca1720489b4f940a`. Reference: `tools/charmm_reference_energy.py`, charmm36_2024.xml, 5TER.

| engine | printed | on CODATA constant |
|---|---|---|
| MOSAICS | 3775.9940254863 | 3775.9959246181 |
| OpenMM | 3775.9959267089 | 3775.9959267089 |

MOSAICS - OpenMM = -2.091e-06 kcal/mol. sander and GROMACS have no CHARMM36 nucleic topology here.
