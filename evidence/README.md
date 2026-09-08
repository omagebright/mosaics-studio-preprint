# Force-field comparison

MOSAICS energies checked against three independent programs, one structure at
a time, with each program's own output file committed next to the input that
produced it. Every number below can be read out of one of those files.

## Result

| force field | structure | largest difference | which term | result |
|---|---|---|---|---|
| OL21/OL3 | 1efs | 1.8e-04 | angle | pass |
| OL15/OL3 | 1efs | 1.8e-04 | angle | pass |
| parmbsc1 | 1efs | 1.8e-04 | angle | pass |
| parmbsc0 | 1efs | 1.8e-04 | angle | pass |
| ff14SB | peptide | 9.6e-06 | angle | pass |
| ff19SB | peptide | 9.6e-06 | angle | pass, see F-06 for the blank GROMACS cell |

In kcal/mol, MOSAICS minus OpenMM, after the electrostatic correction
described in engines.txt. The energies themselves, term by term and program
by program, are in results.txt. The angle is the largest difference on every
row, and the reason is known: it is F-12, and it is an error in the reference
topology, not in MOSAICS.

Before that correction the largest difference on 1efs is electrostatic, and it
is 1.5e-03 kcal/mol: three constants compiled into MOSAICS are rounded. That
is F-13. results.txt prints both values, as each program reports them and
corrected to one constant.

## The five files to read

    results.txt    every energy, term by term, for all four programs, with
                   the file each number came from
    coverage.txt   which force fields and structures are done and which are
                   not, with a reason for every gap
    systems.txt    the structures, what is in them, and the command that
                   builds each one
    engines.txt    the four programs, their versions, how each is run, and
                   the electrostatic constant each one uses
    findings.txt   every problem found, fixed or open, each with an
                   identifier that does not change

## Reading a number out of a file

    sander     awk '/FINAL RESULTS/,0' <structure>/sander/*_sander_sp.out
    OpenMM     cat <structure>/openmm/openmm_sp.out
    GROMACS    tail -1 <structure>/gromacs/sp.xvg
    MOSAICS    grep 'Initial .* energy' <structure>/mosaics/mosaics_sp.out

Two of those files check themselves, and a reader should confirm rather than
take them on trust.

The OpenMM file prints its two halves of the non-bonded energy, their sum, and
the combined value it split. The last two must agree.

The MOSAICS file prints each part of a term next to the total for that term.
Onefour-Coulomb plus Onefour-Lennard-Jones must equal the Onefour total, and
the same holds for Inter and for Torsion. If they do not add up, the structure
moved during the run and the numbers describe two different shapes.

## What counts as passing

Two separate questions are asked, and they are not the same question.

**Does MOSAICS reproduce the force field?** The difference must be no larger
than the reference program can resolve, plus one part in a million of the
number itself. A term further away than that is recorded as unexplained. It is
not called a failure, because the cause may lie in the reference, as F-12
shows. It is never called a pass either, and it gets an identifier in
findings.txt.

**Would the difference change a result?** Measured against 0.002 kcal/mol,
the level Peter has said matters, and against thermal energy at room
temperature, 0.59 kcal/mol. This column never closes an entry. A difference
can be unexplained and inconsequential at the same time. That is a real state
and it is worth saying rather than hiding behind an adjective.

A comparison counts only if all of the following hold.

1. Nothing moved. Every proposal width is zero, and where the run is driven in
   Cartesian coordinates the time step must be zero as well, or the structure
   drifts even with the widths at zero. Do not assume this: check that the
   parts of each term still add up to that term's total.
2. MOSAICS and the reference program read the same coordinates, digit for
   digit. This is not a formality. See F-02.
3. Electrostatic terms are compared only after correcting to one constant, and
   each program's own constant is stated.
4. The reference topology can be rebuilt by a command committed here.
5. Every result names the program it was measured against. Agreeing with one
   program is not a pass.

## Words used here

    term             one part of the energy, such as bond or angle
    topology         the file saying which atoms are bonded and what
                     parameters each carries
    structure file   the file saying where the atoms are
    single point     the energy of one fixed shape, nothing moving
    Onefour          the interaction between two atoms three bonds apart,
                     which every force field treats separately from the rest
    Inter            the interaction between atoms further apart than that
    correction map   a table adjusting a residue's backbone energy according
                     to its two backbone angles. ff19SB and CHARMM36 both use
                     them, and they are why those two are harder to convert
    deck             the set of files holding one force field's parameters in
                     the form MOSAICS reads
