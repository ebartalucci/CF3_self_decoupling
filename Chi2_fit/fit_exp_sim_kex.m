% Least square fit to extract K_ex and T_2
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 19.04.24
% Project: CF3 self decoupling
clear all;
tic;

td=3072;
name= 'spectra/TLA_S_14khz_exp_200_161023/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_14khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_14khz=spectrum_cf3_tla_s_14khz(1:2:end)+1i*spectrum_cf3_tla_s_14khz(2:2:end);

name= 'spectra/TLA_S_30khz_exp_14_121023/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_30khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_30khz=spectrum_cf3_tla_s_30khz(1:2:end)+1i*spectrum_cf3_tla_s_30khz(2:2:end);

name= 'spectra/TLA_S_60khz_exp_15_121023/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_60khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_60khz=spectrum_cf3_tla_s_60khz(1:2:end)+1i*spectrum_cf3_tla_s_60khz(2:2:end);

name= 'spectra/TLA_rac_14khz_exp_10_221123/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_14khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_14khz=spectrum_cf3_tla_rac_14khz(1:2:end)+1i*spectrum_cf3_tla_rac_14khz(2:2:end);

name= 'spectra/TLA_rac_30khz_exp_10_211123/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_30khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_30khz=spectrum_cf3_tla_rac_30khz(1:2:end)+1i*spectrum_cf3_tla_rac_30khz(2:2:end);

name= 'spectra/TLA_rac_60khz_exp_13_211123/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_60khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_60khz=spectrum_cf3_tla_rac_60khz(1:2:end)+1i*spectrum_cf3_tla_rac_60khz(2:2:end);

data = zeros(6,td/2);
data(1,:)=spectrum_cf3_tla_s_14khz;
data(2,:)=spectrum_cf3_tla_s_30khz;
data(3,:)=spectrum_cf3_tla_s_60khz;
data(4,:)=spectrum_cf3_tla_rac_14khz;
data(5,:)=spectrum_cf3_tla_rac_30khz;
data(6,:)=spectrum_cf3_tla_rac_60khz;

datap=zeros(6,32768);
datap(1,:)=phase1(data(1,:),100000,32768,20,66,-38,0,2,15900,67);
datap(2,:)=phase1(data(2,:),100000,32768,20,210,-38,0,2,15900,67);
datap(3,:)=phase1(data(3,:),100000,32768,20,210,-38,0,2,15900,67);
datap(4,:)=phase1(data(4,:),100000,32768,20,62,-38,0,2,15900,67);
datap(5,:)=phase1(data(5,:),100000,32768,20,210,-38,0,2,15900,67);
datap(6,:)=phase1(data(6,:),100000,32768,20,210,-38,0,2,15900,67);

xax=((0:32767)/32768-0.5)*100;
plot(xax,real(datap))


range = 16385:19334;
datapx= datap(:,range);
datapx(1,:)=datap(1,range+819);
xaxp=xax(range);
plot(xaxp,datapx)

sw=xaxp(end)-xaxp(1);
np = length(xaxp);
xax1=((0:np-1)/np-0.5)*sw;
time = (0:np-1)/(sw*1000);

datapx1=real(datapx);
for k=1:6
  datapx1(k,:) = datapx1(k,:)-mean(datapx1(k,1:round(np/4)));
  datapx1(k,:) = datapx1(k,:)/max(datapx1(k,:));
end

plot(xax1,datapx1)

p0 = [-250 400 2 0];
p =zeros(6,4);

options=optimset('MaxFunEvals',10000);
options=optimset(options,'MaxIter',10000);
for k=1:6
  if k>4
    p0(1)=0;
  end
  [p(k,:) resnorm] = lsqcurvefit(@fit_fun, p0, time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
  [p(k,:) resnorm] = lsqcurvefit(@fit_fun, p(k,:), time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
end

datasx1=zeros(size(datapx1));

for k=1:6
  datasx1(k,:) = fit_fun(p(k,:),time);
end

for k=1:6
  subplot(2,3,k)
  plot(xax1,datapx1(k,:),xax1,datasx1(k,:))
  xlabel('\nu [kHz]')
  ylabel('intensity')
  legend('exp','fit')
  line = sprintf('k_{ex} = %5.1f s^{-1}',p(k,2));
  text(1,0.6,line)
  axis([-5 5 -0.1 1.1])
  switch k
  case 1
    title('S 14 kHz')
  case 2
    title('S 30 kHz')
  case 3
    title('S 60 kHz')
  case 4
    title('rac 14 kHz')
  case 5
    title('rac 30 kHz')
  case 6
    title('rac 60 kHz')
  end
end
orient('landscape')

print -dpdf -fillpage figure_fit_cos2_apod.pdf

