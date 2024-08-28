clear all
close all

% Measured MAS frequencies and fields
nurexp=[0];   % in kHz
B0=[500]       % in MHz
nsites=1;       % number of site specific data/resonances

% Preallocate according to number of vmas x B0 conditions
fit2=cell(1);
X2=cell(1);
Y2=cell(1);

%% T1 @ 1.2GHZ
count=0
for n=1:nsites   % loop over no. of resonances
for m=1:length(nurexp)  % loop over MAS or field in case multiple fields/MAS have been measured
    
    count=count+1;

% Extraction of raw data from Topspin output file
file=sprintf('static_t2_rac_TLA_new.txt',nurexp(m))

a=load(file)
Ipeaks=find(a(:,1)==1)

if n==nsites
xdata=a(Ipeaks(n):end,2);
ydata=a(Ipeaks(n):end,3);
else
xdata=a(Ipeaks(1):Ipeaks(1+1)-1,2);
ydata=a(Ipeaks(n):Ipeaks(n+1)-1,3);   
end

% Initial parameters for non-linear fit (Intensity, T1) Set a reasonable 
% value for T1. T1 should be close to our experimental values for good fitting
x0=[ydata(1), 5e-5];   % Rmk: x0(2) start parameter for T1 relaxation time (in s)

       % Functional form of T1 build-up (change in case of T2/T1rho fit into
    fun = @(x,xdata)x(1)*exp(-1/(1*x(2))*xdata)
    %fun = @(x,xdata)x(1)*(1-exp(-1/(1*x(2))*xdata));
 
    xfit = lsqcurvefit(fun,x0,xdata,ydata)
    fit=fun(xfit,xdata);
    T1(n,m)=xfit(2)
      
    % Interpolation for graphical representation of fit function
    xdata_fit2=linspace(0,max(xdata)+1,100);
    fit2{n,m}=fun(xfit,xdata_fit2);
    X2{n,m}=xdata;
    Y2{n,m}=ydata;
    
    % Error calculation via bootstrapping
    resid = ydata-fit;
    foo = fitdist(resid,'normal');
  
   error(n,m,:) = std(bootstrp(200, @(bootr)lsqcurvefit(fun,x0,xdata, fit+bootr),resid))  
    
   % Graphical representation of results (fit vs. raw data). Change number
   % of subplots(m,q,n) according to number of fitted resonances
   fig20=figure(20)
   fig20.Position=[92 72 1280 714];
   if n==1
   plot(xdata,ydata/max(ydata),'ro')
   hold on
   plot(xdata_fit2,fit2{n,m}/max(fit2{n,m}),'k'),
   xlim([0 max(xdata)+1])
   ylim([0 1])
   else
        subplot(5,2,n)
   plot(xdata,ydata/max(ydata),'o')
   hold on
   plot(xdata_fit2,fit2{n,m}/max(fit2{n,m}),'--k'),
   xlim([0 max(xdata)+1])
   ylim([0 1]) 
   end
   
end
end
