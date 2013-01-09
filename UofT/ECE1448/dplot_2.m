
	clear all
	close all
 	clc
	k=input('what data?','s')
    data=load(k);
	idata=data(:,2);
	odata=data(:,3);
    mfp=mean(data(:,4))
    %mfp=4.0E-7
    idata = idata(idata ~= 0);
    odata = odata(odata ~= 0);
    %idata=idata./(mfp);
    %odata=odata./(mfp);
    %[n,xout]=hist(idata,100);
    %[mu,sigma,muci,sigmaci] = normfit(idata)
    %bar(xout, n./sum(n), 1)
    %histnorm(ydata,100)
    %normplot(ydata)
    %histfit(ydata,100,'lognormal')
    [fi,xi]=ksdensity(idata,'support','positive');
    mean(idata)
    plot(xi,fi,':k')
    hold on
    [fo,xo]=ksdensity(odata,'support','positive');
    mean(odata)
    plot(xo,fo,'--k')
    hold on
    x=0:.01*mfp:2*mfp;
    y1=4.*((0.5*mfp)^-2).*x.*exp(-2.*x./(0.5*mfp));
    y2=((0.5*mfp)^-1).*exp(-x./(0.5*mfp));
    plot(x,y1,'k',x,y2,'-.k')
    legend('incoming','outgoing','Goodman-Wachman','MFP distribution')
    xlabel('r-a [m]')
    ylabel('Frequency')

