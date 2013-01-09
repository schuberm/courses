clear all;
clc;


al=dlmread('alfc.txt');
si=dlmread('sifc.txt');

ral=sqrt(al(:,1).^2+al(:,2).^2+al(:,3).^2);
rsi=sqrt(si(:,1).^2+si(:,2).^2+si(:,3).^2);

figure(1);


p1=semilogy(ral(:,1),abs(al(:,4)),'.','MarkerSize',15,'Color','b');hold on;
p2=semilogy(rsi(:,1),abs(si(:,4)),'.','MarkerSize',15,'Color','r');hold on;


legend([p1 p2],'Aluminum','Silicon');xlabel('Distance [lattice constant]'); ylabel('Force Constants [Ry/Bohr^2]');
print('-dpsc','fc.eps');

