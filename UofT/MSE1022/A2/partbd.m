figure(1);
 data=load('conjgrad_steep.dat');
 data1=load('conjgrad_conjg.dat');
 plot(data(:,1),log10(data(:,2)),data1(:,1),log10(data1(:,2))); grid on; xlabel('Number of Iterations');ylabel('Largest Force [eV/Angstrom] (log scale)');
print('-dpsc','partbd.eps');
