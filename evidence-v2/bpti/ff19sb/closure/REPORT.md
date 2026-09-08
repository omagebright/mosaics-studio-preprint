# BPTI ff19SB dihedral term: MOSAICS-current vs OpenMM, closed

System: `paper/evidence-v2/bpti/ff19sb/topology/bpti_ff19sb.prmtop` + `bpti_ff19sb_matched.rst7`
(byte-identical to the copies under `openmm/`). MOSAICS output: `mosaics-current/mosaics_sp.out`.
RTF: `work/pr-fix-staging-nme-rna/params/top_openmm-ff19sb_protein.rtf`.
Deck: `work/pr-fix-staging-nme-rna/params/db_ff19sb/mosaics_openmm-ff19sb.tors_and_impr`.

## Headline

The "remaining 0.0119 kcal/mol" does not exist in the engine. It was a bookkeeping error in the
earlier analysis: residue 58 was resolved against the RTF's `ALA` record instead of `CALA`, which
MOSAICS actually uses (`mosaics_sp.out` line 222: `58. residue is 58  CALA`). That single error
produced both "OXT improper missing" (it is present in `CALA`) and a wrong ALA58 N-improper order.
After correcting it, the whole -0.042080 kcal/mol dihedral gap closes to -8.5e-08.

Reference values (all recomputed by `/tmp/bpti19/closure.py`):

| quantity | proper | improper | sum |
|---|---|---|---|
| OpenMM PeriodicTorsionForce (Reference platform) | 224.973308327 | 5.550584363 | 230.523892690 |
| analytic from prmtop, prmtop phase 3.141594 rad | 224.973308327 | 5.550584363 | 230.523892690 (diff -5.7e-13) |
| MOSAICS printed | 224.973266426 | 5.508546613 | 230.481813039 |
| MOSAICS emulated from RTF + deck (this script) | 224.973266426 (diff -1.7e-10) | 5.508546613 (diff +6.5e-14) | |

## Contributions

| # | contribution | kcal/mol | how computed |
|---|---|---|---|
| I1 | **ASP3 backbone-N improper evaluated in RTF order `-C H N CA` instead of prmtop order `C CA N H`**. Same k = 1.1, same atoms; phi = +158.395 deg (RTF) vs -156.840 deg (prmtop). E = 0.298255178 vs 0.340296072. | **-0.042040894** | `closure.py`: per-improper (built - prmtop) on the 1 improper whose atom order differs. It is the only one: 197 records are the same order as tleap's, and ALA58 uses `CALA`'s `-C CA N H`, which *is* the prmtop order. |
| I2 | Phase: deck encodes exactly pi (p0 = -p2 = 10567.68/1107.09/1006.45), prmtop stores 3.14159400 rad (`%FLAG DIHEDRAL_PHASE`, line 882). 198 other impropers. | +0.000003144 | sum over the 198 same-order impropers of E(pi) - E(3.141594); cross-checked from prmtop k with exact pi: +0.000003144. |
| P1 | Same phase difference over the 3791 proper terms (2391 unique quads). | -0.000041817 | analytic prmtop sum with exact pi/0 phase minus prmtop-phase sum. |
| P2 | Deck coefficient rounding on proper terms (largest single quad 8.9e-10, `H-N3-CX-C8`). | -0.000000085 | `proper_emulate.py`: deck-emulated MOSAICS proper (224.973266426, matches printed to 1.7e-10) minus exact-pi analytic (224.973266510). |
| | **sum I1 + I2 + P1 + P2** | **-0.042079652** | |
| | MOSAICS - OpenMM, dihedral (`energy_table.txt` line 15: -4.2080e-02) | -0.042079652 | |

Deck-k audit (hypothesis 1): 0 of 199 impropers get a deck k different from the prmtop k. Every
RTF-order type quadruple resolves to an entry with the prmtop's k (10.5 / 1.1 / 1.0); no wildcard
shadowing. Improper count: prmtop 199, MOSAICS-built 199, 0 prmtop-only, 0 built-only.
CYX (hypothesis 3): 6 CYX residues carry ordinary `-C H N CA` / `CA +N C O` records, all same
order as prmtop; no disulfide improper exists in either engine. N-terminus: NARG has no `-C` record
and the prmtop has no N3-centred improper on ARG1, so nothing is dropped. C-terminus: CALA's `+N`
record cannot form and the prmtop has no such term either; CALA's `CA O C OXT` is built (0.059813).

## Residual after all identified pieces

- improper: MOSAICS - OpenMM = -0.042037750; I1 + I2 = -0.042037750; **residual -6.7e-14**
- proper: MOSAICS - OpenMM = -0.000041902; P1 + P2 = -0.000041902; **residual < 1e-12** (P1 alone leaves -8.5e-08, which P2 accounts for)
- dihedral total: **residual -8.5e-08 with I1+I2+P1, 0 to printed precision with P2 included**

## Why the earlier analysis had 0.011867 left over

| item | earlier | correct | error |
|---|---|---|---|
| ALA58 `CA-O-C-OXT` | "no RTF record", 0.059815 lost | in `CALA`, built, 0.059813 | -0.059813 |
| ALA58 N improper | ALA order `-C H N CA` -> 0.292025 | CALA order `-C CA N H` -> 0.220346 (= prmtop) | +0.071679 |
| net | | | **+0.011866** |

Earlier "residual" 5.520414 - 5.508546613 = 0.011867. Match to 1.5e-06 (the rest is the 6-decimal
rounding in the 5.520414 figure). Nothing remains.

## Reproduce

```
E=/Users/bright/Documents/MOSAICS/paper/evidence-v2/bpti/ff19sb
P=/Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/params
PY=/opt/homebrew/Caskroom/miniforge/base/envs/mosaics-ff-local/bin/python
cd /tmp/bpti19
$PY closure.py        $E/topology/bpti_ff19sb.prmtop $E/topology/bpti_ff19sb_matched.rst7 $P/top_openmm-ff19sb_protein.rtf $P/db_ff19sb/mosaics_openmm-ff19sb.tors_and_impr $E/mosaics-current/mosaics_sp.out
$PY proper_emulate.py $E/topology/bpti_ff19sb.prmtop $E/topology/bpti_ff19sb_matched.rst7 $P/top_openmm-ff19sb_protein.rtf $P/db_ff19sb/mosaics_openmm-ff19sb.tors_and_impr
```

Outputs saved: `/tmp/bpti19/closure.txt`, `/tmp/bpti19/proper.txt`. Supporting scripts:
`locate_residual.py` (per-improper table, 199 rows), `perterm_openmm.py` (per-term OpenMM vs
analytic, 3990 rows, max |dE| 6.6e-06 all from the 3.141594 phase), `openmm_reference.py`,
`signcheck.py` (confirms the dihedral sign convention matches OpenMM on a 4-atom test).

## What was and was not verified

- Verified: OpenMM 8.5.2 Reference platform reproduces the analytic prmtop sum to 5.7e-13; the
  RTF+deck emulation of MOSAICS (tors.cpp:186-210 lookup, tors.cpp:756-779 energy, read_rtf.cpp
  terminal-residue naming) reproduces the printed MOSAICS improper to 6.5e-14 and proper to 1.7e-10.
- Not done: MOSAICS was not re-run; the engine's own per-improper list was not extracted. The
  claim that MOSAICS builds ASP3's improper as `-C H N CA` rests on the RTF record plus the
  6.5e-14 agreement of the emulation with the printed total, not on an engine probe.
- No file under `~/Documents/MOSAICS/engine` or `~/Documents/MOSAICS/paper` was modified.
