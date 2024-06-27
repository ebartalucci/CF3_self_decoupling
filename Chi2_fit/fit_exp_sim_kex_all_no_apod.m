% Least square fit to extract K_ex and T_2
% Include all 14 spectra with increasing MAS
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Use cos^2 apodization function in processing 
% Last update: Aachen, 17.06.24
% Project: CF3 self decoupling

clear all;
tic;

% Printing file
fileID = fopen('results_fit_exp_sim_kex_all_no_apod.txt','w');

%% Load 14 spectra

td=3072;

% S-TFLA
name= 'spectra/13C_CP_MAS_dependent_TLA_S_14khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_14khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_14khz=spectrum_cf3_tla_s_14khz(1:2:end)+1i*spectrum_cf3_tla_s_14khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_17p5khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_17p5khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_17p5khz=spectrum_cf3_tla_s_17p5khz(1:2:end)+1i*spectrum_cf3_tla_s_17p5khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_22khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_22khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_22khz=spectrum_cf3_tla_s_22khz(1:2:end)+1i*spectrum_cf3_tla_s_22khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_30khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_30khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_30khz=spectrum_cf3_tla_s_30khz(1:2:end)+1i*spectrum_cf3_tla_s_30khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_40khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_40khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_40khz=spectrum_cf3_tla_s_40khz(1:2:end)+1i*spectrum_cf3_tla_s_40khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_50khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_50khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_50khz=spectrum_cf3_tla_s_50khz(1:2:end)+1i*spectrum_cf3_tla_s_50khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_S_60khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_60khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_s_60khz=spectrum_cf3_tla_s_60khz(1:2:end)+1i*spectrum_cf3_tla_s_60khz(2:2:end);

% Rac-TFLA
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_14khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_14khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_14khz=spectrum_cf3_tla_rac_14khz(1:2:end)+1i*spectrum_cf3_tla_rac_14khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_17p5khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_17p5khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_17p5khz=spectrum_cf3_tla_rac_17p5khz(1:2:end)+1i*spectrum_cf3_tla_rac_17p5khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_22khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_22khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_22khz=spectrum_cf3_tla_rac_22khz(1:2:end)+1i*spectrum_cf3_tla_rac_22khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_30khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_30khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_30khz=spectrum_cf3_tla_rac_30khz(1:2:end)+1i*spectrum_cf3_tla_rac_30khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_40khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_40khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_40khz=spectrum_cf3_tla_rac_40khz(1:2:end)+1i*spectrum_cf3_tla_rac_40khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_50khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_50khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_50khz=spectrum_cf3_tla_rac_50khz(1:2:end)+1i*spectrum_cf3_tla_rac_50khz(2:2:end);

name= 'spectra/13C_CP_MAS_dependent_TLA_rac_60khz/fid';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_60khz = fread(fid,td,'int');
fclose(fid);
spectrum_cf3_tla_rac_60khz=spectrum_cf3_tla_rac_60khz(1:2:end)+1i*spectrum_cf3_tla_rac_60khz(2:2:end);

% Procs
data = zeros(14,td/2);
data(1,:)=spectrum_cf3_tla_s_14khz;
data(2,:)=spectrum_cf3_tla_s_17p5khz;
data(3,:)=spectrum_cf3_tla_s_22khz;
data(4,:)=spectrum_cf3_tla_s_30khz;
data(5,:)=spectrum_cf3_tla_s_40khz;
data(6,:)=spectrum_cf3_tla_s_50khz;
data(7,:)=spectrum_cf3_tla_s_60khz;

data(8,:)=spectrum_cf3_tla_rac_14khz;
data(9,:)=spectrum_cf3_tla_rac_17p5khz;
data(10,:)=spectrum_cf3_tla_rac_22khz;
data(11,:)=spectrum_cf3_tla_rac_30khz;
data(12,:)=spectrum_cf3_tla_rac_40khz;
data(13,:)=spectrum_cf3_tla_rac_50khz;
data(14,:)=spectrum_cf3_tla_rac_60khz;

datap=zeros(14,32768);
% S-TFLA
datap(1,:)=proc_fid(data(1,:),100000,32768,0,66,-38,0,0,15900,67);
datap(2,:)=proc_fid(data(2,:),100000,32768,0,78,-38,0,0,15900,67);
datap(3,:)=proc_fid(data(3,:),100000,32768,0,64,-38,0,0,15900,67);
datap(4,:)=proc_fid(data(4,:),100000,32768,0,210,-38,0,0,15900,67);
datap(5,:)=proc_fid(data(5,:),100000,32768,0,210,-38,0,0,15900,67);
datap(6,:)=proc_fid(data(6,:),100000,32768,0,210,-38,0,0,15900,67);
datap(7,:)=proc_fid(data(7,:),100000,32768,0,210,-38,0,0,15900,67);

% rac-TFLA
datap(8,:)=proc_fid(data(8,:),100000,32768,0,62,-38,0,0,15900,67);
datap(9,:)=proc_fid(data(9,:),100000,32768,0,60,-38,0,0,15900,67);
datap(10,:)=proc_fid(data(10,:),100000,32768,0,60,-38,0,0,15900,67);
datap(11,:)=proc_fid(data(11,:),100000,32768,0,210,-38,0,0,15900,67);
datap(12,:)=proc_fid(data(12,:),100000,32768,0,206,-38,0,0,15900,67);
datap(13,:)=proc_fid(data(13,:),100000,32768,0,210,-38,0,0,15900,67);
datap(14,:)=proc_fid(data(14,:),100000,32768,0,210,-38,0,0,15900,67);

xax=((0:32767)/32768-0.5)*100;
plot(xax,real(datap))

range = 16385:19334;
datapx= datap(:,range);
xaxp=xax(range);

% shift spectra
datapx(1,:)=datap(1,range+819); % shift experiment to match with the others
datapx(2,:)=datap(2,range+819); % shift experiment to match with the others
datapx(3,:)=datap(3,range+819); % shift experiment to match with the others
datapx(10,:)=datap(10,range-30); % shift experiment to match with the others
datapx(11,:)=datap(11,range-30); % shift experiment to match with the others
%datapx(12,:)=datap(12,range-30); % shift experiment to match with the others

plot(xaxp,datapx)

sw=xaxp(end)-xaxp(1);
np = length(xaxp);
xax1=((0:np-1)/np-0.5)*sw;
time = (0:np-1)/(sw*1000);

datapx1=real(datapx);

for k=1:14
  datapx1(k,:) = datapx1(k,:)-mean(datapx1(k,1:round(np/4)));
  datapx1(k,:) = datapx1(k,:)/max(datapx1(k,:));
end

plot(xax1,datapx1)

p0 = [-250 400 2 0];
p =zeros(6,4);

options=optimset('MaxFunEvals',10000);
options=optimset(options,'MaxIter',10000);
for k=1:14
  if k>8
    p0(1)=0;
  end
  [p(k,:) resnorm] = lsqcurvefit(@fit_fun, p0, time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
  [p(k,:) resnorm] = lsqcurvefit(@fit_fun, p(k,:), time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
end

datasx1=zeros(size(datapx1));

for k=1:14
  datasx1(k,:) = fit_fun(p(k,:),time);
end

for k=1:14
    subplot(7,2,k)
    plot(xax1,datapx1(k,:),xax1,datasx1(k,:))
    xlabel('\nu [kHz]')
    ylabel('intensity')
    legend('exp','fit')
    line = sprintf('k_{ex} = %5.1f s^{-1}',p(k,2));

    % Write the results to the text file
    fprintf(fileID, 'Spectrum: %f\n', k);
    fprintf(fileID, 'k_{ex} = %5.1f s^{-1}',p(k,2));

    text(1,0.6,line)
    axis([-2 2 -0.1 1.2])

    switch k
        case 1
            title('S-TFLA 14 kHz')
        case 2
            title('S-TFLA 17.5 kHz')
        case 3
            title('S-TFLA 22 kHz')
        case 4
            title('S-TFLA 30 kHz')
        case 5
            title('S-TFLA 40 kHz')
        case 6
            title('S-TFLA 50 kHz')
        case 7
            title('S-TFLA 60 kHz')
        case 8
            title('rac-TFLA 14 kHz')
        case 9
            title('rac-TFLA 17.5 kHz')
        case 10
            title('rac-TFLA 22 kHz')
        case 11
            title('rac-TFLA 30 kHz')
        case 12
            title('rac-TFLA 40 kHz')
        case 13
            title('rac-TFLA 50 kHz')
        case 14
            title('rac-TFLA 60 kHz')
    end
end

% Close the text file
fclose(fileID);

print -dpdf -fillpage output/stacked_figure_all_fits_no_apod.pdf

