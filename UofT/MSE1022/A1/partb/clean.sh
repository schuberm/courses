#!/bin/bash
grep -A 1 "Step PotEng KinEng TotEng Press Volume" log.lammps > si1.txt
sed -e ''s/[a-zA-Z]//g'' -e ''s/--//g'' -e 's/^[ \t]*//g' -e ''/^$/d'' si1.txt > si2.txt
awk '{ getline ln < "a0.txt"; print ln" "$0 }' si2.txt > si.txt
rm si1.txt
rm si2.txt
