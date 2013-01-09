# gnuplot script to generate plots for thermodynamic output
# Format of the thermo output (from lammps script or default) is:
# Step Atoms Temp PotEng Lx Ly Lz Pzz Press
set term postscript enhanced # eps file format
set output 'temp_plot.eps'
set xlabel 'time (ps) = Step*dt'
set ylabel 'Temperature (K)'
set grid
dt=0.005
plot 'thermo.log' u ($1*dt):3 title 'temperature' w l

