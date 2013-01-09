# gnuplot script for generating stress-strain curve
set term postscript enhanced
set output 'se_curve.eps'
set xlabel 'Strain'
set ylabel 'Stress (MPa)'
set grid 
erate=0.00001
plot 'mainthermo.log' u ($1*erate):(-$7/10) title 'Stress-Strain curve' w l  

