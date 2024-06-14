% Test chisquare fit to extract K_ex and T_2 using matthias 1D fits
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Project: CF3 self decoupling

clear all;
tic;

%% Get spectra and process 
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

% phasing and processing using kathrin's function
datap=zeros(6,32768);
datap(1,:)=proc_fid(data(1,:),100000,32768,20,66,-38,0,2,15900,67);
datap(2,:)=proc_fid(data(2,:),100000,32768,20,210,-38,0,2,15900,67);
datap(3,:)=proc_fid(data(3,:),100000,32768,20,210,-38,0,2,15900,67);
datap(4,:)=proc_fid(data(4,:),100000,32768,20,62,-38,0,2,15900,67);
datap(5,:)=proc_fid(data(5,:),100000,32768,20,210,-38,0,2,15900,67);
datap(6,:)=proc_fid(data(6,:),100000,32768,20,210,-38,0,2,15900,67);

% test to see if things fit
xax=((0:32767)/32768-0.5)*100;
plot(xax,real(datap))

% i think this is the zooming in CF3 region
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

% phases
p0 = [-250 400 2 0];
p =zeros(6,4);

options=optimset('MaxFunEvals',10000);
options=optimset(options,'MaxIter',10000);

for k=1:6
  if k>4
    p0(1)=0;
  end
  [p(k,:), resnorm] = lsqcurvefit(@fit_fun, p0, time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
  [p(k,:), resnorm] = lsqcurvefit(@fit_fun, p(k,:), time,datapx1(k,:),[-10000 0 0 -10000],[10000 10000 0.1 10000],options);
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
    title('({\itS})-TFLA 14kHz')
  case 2
    title('({\itS})-TFLA 30kHz')
  case 3
    title('({\itS})-TFLA 60kHz')
  case 4
    title('({\itrac})-TFLA 14kHz')
  case 5
    title('({\itrac})-TFLA 30kHz')
  case 6
    title('({\itrac})-TFLA 60kHz')
  end
end
orient('landscape')

print -dpdf -fillpage figure_fit_cos2_apod.pdf


%% Chi square simulations
J_cf = 280; 
fileID = fopen('chisquare_mins.txt','w');


% Define different t values that fits the experimental sizes
t_values = {linspace(0, 0.1, 8192), linspace(0, 0.1, 8192), ...
            linspace(0, 0.1, 131072), linspace(0, 0.1, 262144), ...
            linspace(0, 0.1, 8192), linspace(0, 0.1, 131072)};

% loop on exp spectra
for k=1:6
    switch k
        case 1 % s 14khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(500, 1000, 50);
            T_2_values = linspace(0.01, 0.1, 50); 
            offset = 35600;
        
        case 2 % s 30khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(1, 300, 50);
            T_2_values = linspace(0.01, 0.1, 50);
            offset = 22300;

        case 3 % s 60 khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(1, 100, 50);
            T_2_values = linspace(0.01, 0.1, 50);
            offset = 1;

        case 4 % rac 14khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(1, 500, 50);
            T_2_values = linspace(0.01, 0.1, 50);
            offset = 1;

        case 5 % rac 30 khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(1, 200, 50);
            T_2_values = linspace(0.01, 0.1, 50);
            offset = 24200;

        case 6 % rac 60 khz
            current_spectrum = datapx1(k,:);
            k_ex_values = linspace(1, 100, 50);
            T_2_values = linspace(0.01, 0.1, 50);
            offset = 1;
    end
    
    % NMR variables
    t = t_values{t_index}; % time domain size (s) need to match experiment
    DW = t(2) - t(1); % dwell time
    SW = 1/(2*DW); % spectral width (Hz)
    NP = length(zeros(size(t))); % Length of the signal
      
    % Store chi-square statistics
    chi_square_matrix = zeros(length(k_ex_values), length(T_2_values));

    % Compute the  frequency axis
    f = linspace(-SW, SW, NP); % Frequency axis
    f2 = linspace(-SW+offset, SW+offset, NP); % Shifted axis by offset
    
    % Compute final spectrum for all combinations of k_ex and T_2 values
    for j = 1:length(k_ex_values)
        for T = 1:length(T_2_values)
            
            k_ex = k_ex_values(j);
            T_2 = T_2_values(T);

            % Modified Bloch equtions
            L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
                 k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit
            
            L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
                  3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
                  0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
                  0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];
        
            % Compute the matrix exponentials to get the time evolutions
            U1 = zeros(2, 2, length(time)); % Preallocate U1 matrix
            for i = 1:length(time)
                U1(:,:,i) = expm(L1 * abs(time(i)));
            end
        
            U2 = zeros(4, 4, length(time)); % Preallocate U2 matrix
            for i = 1:length(time)
                U2(:,:,i) = expm(L2 * abs(time(i)));
            end
        
            % Get the signals
            vec_u1 = [2; 2];
            vec_u2 = ones(4, 1);
        
            % Compute signal 1
            signal1 = zeros(size(time));
            for i = 1:length(time)
                signal1(i) = sum(sum(U1(:,:,i) .* vec_u1));
            end
            
            % Compute signal2
            signal2 = zeros(size(time));
            for i = 1:length(time)
                signal2(i) = sum(sum(U2(:,:,i) .* vec_u2));
            end
        
            % Fourier transform and shift 
            spectrum1 = fft(signal1);
            spectrum1 = fftshift(spectrum1);
            
            spectrum2 = fft(signal2);
            spectrum2 = fftshift(spectrum2);
        
            % Compute final spectrum
            final_spectrum = spectrum1 + spectrum2;

            final_spectrum = real(final_spectrum(:))/max(real(final_spectrum));
            plot(f,final_spectrum,f2,current_spectrum)
            [a1,b1a] = min(abs(f+3000));
            [a1,b1b] = min(abs(f-3000));
            [a1,b2a] = min(abs(f2+3000));
            b2b=b2a+b1b-b1a;
    		plot(f(b1a:b1b),final_spectrum(b1a:b1b),f2(b2a:b2b),current_spectrum(b2a:b2b))
            chi_square_matrix(k, T) = sum((final_spectrum(b1a:b1b)-current_spectrum(b2a:b2b)).^2);
        end
    end

    % Plot chi-square statistics
    % Compute the indices of the minimum chi-square value
    [min_chi_square, min_index] = min(chi_square_matrix(:));
    [min_k_ex_index, min_T_2_index] = ind2sub(size(chi_square_matrix), min_index);
    
    % Retrieve the corresponding k_ex and T_2 values
    min_k_ex = k_ex_values(min_k_ex_index);
    min_T_2 = T_2_values(min_T_2_index);
    
    % Print the results to the console
    disp(['Minimal Chi-square for TLA Spectrum ', num2str(k), ':']);
    disp(['T2 Value: ', num2str(min_T_2)]);
    disp(['k_ex Value: ', num2str(min_k_ex)]);
    disp(['Chi-square Value: ', num2str(min_chi_square)]);
    disp(' ');
    
    % Write the results to the text file
    fprintf(fileID, 'Minimal Chi-square for TLA Spectrum %d:\n', k);
    fprintf(fileID, 'T2 Value: %f\n', min_T_2);
    fprintf(fileID, 'k_ex Value: %f\n', min_k_ex);
    fprintf(fileID, 'Chi-square Value: %f\n\n', min_chi_square);

    % Normalize the chi-square matrix to the range [0, 1]
    normalized_chi_square_matrix = chi_square_matrix / max(max(chi_square_matrix));

    % Plot the chi-square statistics as a heatmap with a specified colormap
    figure(1000);
    subplot(2, 3, k); 
    imagesc(T_2_values, k_ex_values, normalized_chi_square_matrix);
    colormap(jet); 
    colorbar; 
    xlabel('T2');
    ylabel('k_ex');

end



