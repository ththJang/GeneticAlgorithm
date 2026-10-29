#!/bin/bash
# Taehwan Jang jthgreene@kaist.ac.kr
# Genetic-algorithm fit of force-field parameters against reference energies.
# Run from an example directory (e.g. examples/water_rexpon):
#   LMP=/path/to/lmp_mpi /path/to/scripts/run_ga.sh [NGEN]
#
# required files in the current directory
#   ini_chromosomes : initial population, one chromosome per line (population % 4 == 0);
#                     the mutation range is set by mutation_range in min_genetic2.py
#   base.ffield.fit : force-field template with placeholders (see write_ffield.sh)
#   control         : RexPoN control file
#   reffit          : reference energies
#   clusters.dat    : cluster names and subtracted reference energies
#   template/       : LAMMPS inputs and data files copied into every work directory
# output
#   chromosomes     : current population (overwritten every generation)
#   min_history     : best chromosome + its error, one line per generation

SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NGEN=${1:-40000}
NPOP=`grep -c . ini_chromosomes`

echo "min" >> min_history
cp ini_chromosomes chromosomes
for i in `seq 1 $NGEN`
do
	for j in `seq 1 $NPOP`;do
		(cp -r template seq.${j}/; awk -v nr=$j '{if(NR==nr) print $0}' chromosomes > seq.${j}/chromosome;) &
	done
	wait
	for j in `seq 1 $NPOP`;do
		(cd seq.${j}; $SCRIPTS/write_ffield.sh chromosome ../base.ffield.fit test.ffield.RexPoN; ln -s ../control .; $SCRIPTS/evaluate_chromosome.sh ../reffit ../clusters.dat; cd ..;) &
	done
	wait
	rm -f err
	for j in `seq 1 $NPOP`
	do
		cat seq.${j}/err >> err
	done
	rm -r seq.*
	paste chromosomes err > chromosomes2
	python3 $SCRIPTS/min_genetic2.py chromosomes2 > tempo
	head -n1 tempo >> min_history
	tail -n $NPOP tempo > chromosomes
done
