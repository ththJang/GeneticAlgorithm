# GeneticAlgorithm
A framework for a genetic algorithm can be adapted to optimize sequences of numbers.
Written by Taehwan Jang (jthgreene@kaist.ac.kr)

Example:
# Genetic-algorithm fitting of force-field parameters

Fits a set of force-field parameters (a "chromosome") with a genetic algorithm (GA).
The fit minimizes the mean absolute error (MAE) between energies computed with LAMMPS
and reference energies. The example fits 12 lone-pair / hydrogen-bond parameters of
a RexPoN water force field (with PQEq charges) to two kinds of reference data:
water-dimer potential-energy scans and binding energies of water clusters (n = 2–19).


## Layout

```
scripts/
  run_ga.sh                GA driver: evaluates the population and produces the next generation
  min_genetic2.py          selection, crossover and mutation in one generation
  write_ffield.sh          puts one chromosome into the force-field template
  evaluate_chromosome.sh   runs LAMMPS for one force field and writes its MAE to ./err
examples/water_rexpon/
  ini_chromosomes          initial population (12 identical rows)
  base.ffield.fit          RexPoN force field with placeholders lp1..lp4, hb1..hb8
  control                  RexPoN control file
  reffit                   172 reference energies (kcal/mol), see "Reference data"
  clusters.dat             cluster names and the energy subtracted from each
  template/                LAMMPS inputs and data files, copied into each work directory
  slurm.sh                 example batch script
```

## Requirements

- bash, awk, sed, Python ≥ 3.6 (standard library only)
- A LAMMPS build that includes `pair_style reaxff/rexpon`, `pair_style coul/pqeqgauss`,
  `fix psq` and `atom_style pqeq` (RexPoN/PQEq build, LAMMPS 24Jun2022 was used).
  Pass it in the `LMP` environment variable.

## Running

```bash
cd examples/water_rexpon
LMP=/path/to/lmp_mpi ../../scripts/run_ga.sh 40000     # number of generations
```

or submit `slurm.sh` after setting `LMP` and the partition. All chromosomes of a
generation run in parallel as background processes. Each one gets its own
`seq.N/` directory, which is deleted at the end of the generation.

Output in the run directory:

| file          | content                                                             |
|---------------|---------------------------------------------------------------------|
| `min_history` | one line per generation: best chromosome followed by its MAE        |
| `chromosomes` | current population (restart a run by copying it to `ini_chromosomes`) |
| `chromosomes2`| population of the last generation with the MAE as the last column   |

The best parameter set is the line of `min_history` with the smallest last column:
`grep -v '^min' min_history | sort -g -k13 | head -1`.

## Input formats

**Chromosome**: one line of whitespace-separated numbers. The i-th value replaces
the i-th placeholder listed in `PLACEHOLDERS` in `write_ffield.sh`
(`lp1 lp2 lp3 lp4 hb1 ... hb8`). These placeholders appear in `base.ffield.fit`.
To fit other parameters, put new placeholder names in the template and edit
`PLACEHOLDERS`.

**ini_chromosomes**: the initial population, one chromosome per line. The population
size must be a multiple of 4; this is not checked.

**reffit**: one reference energy per line, in this order:

| lines   | data                                                        |
|---------|-------------------------------------------------------------|
| 1–41    | scan1: dimer, monomer 2 translated in 41 steps              |
| 42–68   | scan2: dimer, monomer 2 rotated in 27 steps of 5°           |
| 69–99   | scan3: dimer, monomer 2 translated in 31 steps              |
| 100–136 | scan4: dimer, monomer 2 rotated in 37 steps of 10°          |
| 137–172 | the 36 structures of `clusters.dat`, in the same order      |

The scans are rigid (`run 0` at each step). The computed scan energy is
`E_dimer - MONOMERSUM` (`MONOMERSUM = -541.54` kcal/mol, set in
`evaluate_chromosome.sh`). The computed cluster energy is `E(data.<name>)` minus the
second column of `clusters.dat`. Those values are the FF energies of the separated
monomers: `-270.7742 × n` for the n-mer clusters, and individual values for the
dimers `str01`–`str10`. Both references are fixed numbers. They are **not**
recomputed with the trial parameters.

## Algorithm (`min_genetic2.py`)

`run_ga.sh` calls it as `python3 min_genetic2.py chromosomes2 > tempo`. It needs
Python 3: its first line points to `/usr/bin/python`, which can be Python 2, so
don't run it directly. It prints the best chromosome (with its error) followed by
the new population.

For a population of N chromosomes with errors:

1. Sort by error and keep the best N/2.
2. **Selection**: rank-weighted roulette. Survivor of rank i (0 = best) is drawn with
   weight N/2 − i; the best one gets an extra N/2.
3. **Crossover**: duplicate the selection. Pairs of the copies are blended gene by gene
   with a random β ∈ [0,1): `x1' = β x1 + (1-β) x2`, `x2' = (1-β) x1 + β x2`.
4. **Mutation**: every gene of every chromosome is replaced, with probability 0.1,
   by `ref × U(0, 2)`, where `ref` is that gene in the list `mutation_range`
   (probability: `do_mutate` in `mutation()`).

`mutation_range` is hard-coded in `min_genetic2.py` and equals the first row of
`ini_chromosomes`. If you change the initial values, or fit a different number of
parameters, edit `mutation_range` to match. The script does not check this: with
a different number of genes it silently uses gene `j % len(mutation_range)`'s range.

Consequences worth knowing:
- A gene whose reference value is 0 (here `hb7`) stays 0 for the whole fit.
- A negative reference (here `hb1 = -43.31`) gives a range of [2·ref, 0].
- There is no strict elitism: the best chromosome is usually selected, but it can be
  changed by crossover or mutation, so the best error in `min_history` is not
  guaranteed to decrease monotonically.
- If a LAMMPS evaluation fails or prints a different number of energies than `reffit`
  has, that chromosome gets error 1e10 and is discarded.

