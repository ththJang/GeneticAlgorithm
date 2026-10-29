#!/bin/bash
# Evaluate one force field (test.ffield.RexPoN in the current directory) with LAMMPS
# and write its mean absolute error against the reference energies to ./err
# usage (run inside a work directory holding the template files):
#   evaluate_chromosome.sh REFFIT CLUSTERS
#     REFFIT   : reference energies, one per line, in the order of energyfit1 below
#     CLUSTERS : two columns, <name> <energy subtracted from E(data.<name>)>
# environment: LMP = LAMMPS executable with RexPoN / PQEq (default: lmp_mpi)
#
# energyfit1 = [dimer scans: E(scan i) - MONOMERSUM for each frame of scan1..4.in]
#            + [clusters  : E(data.<name>) - reference column of CLUSTERS]

lammps=${LMP:-lmp_mpi}
reffit=$1
clusters=$2
monomersum=-541.54   # energy of two separated monomers (kcal/mol), subtracted from the scans

rm -f energyfit en_tot
for i in 1 2 3 4; do
	$lammps -in scan${i}.in > scan${i}.log
	grep 'TotEng' scan${i}.log | awk -v mono=$monomersum '{print $9-mono}' > en${i}
	cat en${i} >> energyfit
done
for name in `awk '{print $1}' $clusters`; do
	awk -v name=$name '{if ($1=="read_data") {print $1,"data."name} else {print $0}}' cluster.in > run.in
	$lammps -in run.in > $name.out
	grep TotEng $name.out | awk '{print $9}' >> en_tot
done
awk '{print $2}' $clusters | paste en_tot - > en_tot+ref
awk '{print $1-$2}' en_tot+ref > en_final
cat energyfit en_final > energyfit1

# a failed LAMMPS run would shift the rows; rank such a chromosome last
if [ `wc -l < energyfit1` -ne `wc -l < $reffit` ] || [ `wc -l < en_tot` -ne `wc -l < $clusters` ]; then
	echo "evaluate_chromosome.sh: energy count mismatch in `pwd`" >&2
	echo 1e10 > err
	exit 1
fi
paste $reffit energyfit1 | awk '{d=$1-$2; if(d<0) d=-d; sum+=d} END {print sum/NR}' > err # MAE
