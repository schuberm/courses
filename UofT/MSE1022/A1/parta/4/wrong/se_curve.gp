# gnuplot script for generating stress-strain curve
set term postscript enhanced
set output 'se_curve4.eps'
set xlabel 'Strain'
set ylabel 'Stress (MPa)'
set grid 
erate=0.00001
f1(x)=E1*x
f2(x)=E2*x
#
set style line 1 lt 1 lw 2 pt 3 linecolor rgb "red"
set style line 2 lt 1 lw 2 pt 3 linecolor rgb "green"
set style line 3 lt 1 lw 1 pt 3 linecolor rgb "blue"
#

fit f1(x) '< head -80 eam.txt' u ($1*erate):(-$7/10) via E1
fit f2(x) '< head -80 lj.txt' u ($1*erate):(-$7/10) via E2

#plot '< head -100 eam.txt' u ($1*erate):(-$7/10) title 'Strain Rate EAM' w l ls 1,f1(x) title 'Linear Fit EAM' ls 1,'< head -100 lj.txt' u ($1*erate):(-$7/10) title 'Strain Rate LJ' w l ls 2,f2(x) title 'Linear Fit LJ' ls 2

plot '< head -100 eam.txt' u ($1*erate):(-$7/10) title 'Strain Rate EAM' w l ls 1, '< head -100 lj.txt' u ($1*erate):(-$7/10) title 'Strain Rate LJ' w l ls 2
