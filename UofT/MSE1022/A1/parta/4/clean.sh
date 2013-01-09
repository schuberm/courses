#!/bin/bash
grep -A 1 "Step PotEng KinEng TotEng Press Volume" log.lammps > lj1.txt
sed -e ''s/[a-zA-Z]//g'' -e ''s/--//g'' -e 's/^[ \t]*//g' -e ''/^$/d'' lj1.txt > lj2.txt
awk '{ getline ln < "a0.txt"; print ln" "$0 }' lj2.txt > lj.txt
rm lj1.txt
rm lj2.txt
