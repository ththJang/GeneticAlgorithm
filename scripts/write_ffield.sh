#!/bin/bash
# Substitute one chromosome into the force-field template.
# usage: write_ffield.sh CHROMOSOME_FILE BASE_FFIELD OUTPUT_FFIELD
# The i-th value of the chromosome replaces the i-th placeholder in PLACEHOLDERS.

PLACEHOLDERS=(lp1 lp2 lp3 lp4 hb1 hb2 hb3 hb4 hb5 hb6 hb7 hb8)

chromosome=(`cat $1`)
if [ ${#chromosome[@]} -ne ${#PLACEHOLDERS[@]} ]; then
	echo "write_ffield.sh: chromosome has ${#chromosome[@]} values, expected ${#PLACEHOLDERS[@]}" >&2
	exit 1
fi
cp $2 $3
for ((k=0;k<${#PLACEHOLDERS[@]};k++)); do
	sed -i "s/${PLACEHOLDERS[$k]}/${chromosome[$k]}/" $3
done
