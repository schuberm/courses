# gnuplot script for generating stress-strain curve
set term postscript enhanced
set output 'disp_curve.eps'
set xlabel 'Distance between atoms [A]'
set ylabel 'Potential Energy [eV]'
set grid 
#
set style line 1 lt 1 lw 2 pt 3 linecolor rgb "red"
set style line 2 lt 1 lw 2 pt 3 linecolor rgb "green"
set style line 3 lt 1 lw 2 pt 3 linecolor rgb "blue"

epsilon=0.458 
sigma=2.569

#f(x)=epsilon*((sigma/x)^12-(sigma/x)^6)
#f(x)=4*0.458*((2.569/x)**(12)-(2.569/x)**6)
#set xrange [1:5]
#set yrange [-20:20]

#plot f(x) title 'LJ' ls 1, 'eamclean.txt' u ($1):(-$2-0.5*(27.2*0.529*$3*$3)/($1)) title 'EAM' w l ls 2
#plot 'eamclean.txt' u ($1):($2+0.5*(27.2*0.529*$3*$3)) title 'EAM' w l ls 2
#plot f(x) title 'LJ' ls 1, 'eamclean.txt' u ($1):($2+0.5*(27.2*0.529*$3*$3)) title 'EAM' w l ls 2

plot 'eam.txt' u ($1):($3/108) title 'EAM' w l ls 1, 'lj.txt' u ($1):($3/108) title 'LJ' w l ls 2
