# gnuplot script for generating stress-strain curve
set term postscript enhanced
set output 'pe_curve.eps'
set xlabel 'Lattice Constant [A]'
set ylabel 'Energy [eV]'
set grid 

set style line 1 lt 1 lw 2 pt 3 linecolor rgb "red"
set style line 2 lt 1 lw 2 pt 3 linecolor rgb "green"
set style line 3 lt 1 lw 1 pt 3 linecolor rgb "blue"

f1(x)=a1 + b1*x + c1*x**2
f2(x)=a2 + b2*x + c2*x**2 + d2*x**3 + e2*x**4 + f2*x**5
fit f1(x) 'si.txt' u ($1):($3/216) via a1,b1,c1
fit f2(x) 'si.txt' u ($1):($3/216) via a2,b2,c2,d2,e2,f2

plot 'si.txt' every 100 u ($1):($3/216) title 'Potential Energy' w p, f1(x) title 'Quadratic Fit' ls 1, f2(x) title 'Quintic Fit' ls 2 
 

