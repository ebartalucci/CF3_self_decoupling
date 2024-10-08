% Compare experiments and simulations with best T2 and k_ex
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 18.04.24
% Project: CF3 self decoupling

%% Experimental section
clear all;

% S-TFLA
si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_14khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_14khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_17p5khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_17p5khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_22khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_22khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_30khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_30khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_40khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_40khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_50khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_50khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_S_60khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_60khz = fread(fid,si,'int');
fclose(fid);

% Rac-TFLA
si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_14khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_14khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_17p5khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_17p5khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_22khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_22khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_30khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_30khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_40khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_40khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_50khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_50khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/13C_CP_MAS_dependent_TLA_rac_60khz/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_60khz = fread(fid,si,'int');
fclose(fid);

% Normalize experimental spectral intensity of CF3 peak between 0 and 1
% S-TFLA
spectrum_cf3_tla_s_14khz = spectrum_cf3_tla_s_14khz / max(spectrum_cf3_tla_s_14khz(3500:4000));
spectrum_cf3_tla_s_17p5khz = spectrum_cf3_tla_s_17p5khz / max(spectrum_cf3_tla_s_17p5khz(3500:4000));
spectrum_cf3_tla_s_22khz = spectrum_cf3_tla_s_22khz / max(spectrum_cf3_tla_s_22khz(3500:4000));
spectrum_cf3_tla_s_30khz = spectrum_cf3_tla_s_30khz / max(spectrum_cf3_tla_s_30khz(3500:4000));
spectrum_cf3_tla_s_40khz = spectrum_cf3_tla_s_40khz / max(spectrum_cf3_tla_s_40khz(3500:4000));
spectrum_cf3_tla_s_50khz = spectrum_cf3_tla_s_50khz / max(spectrum_cf3_tla_s_50khz(3500:4000));
spectrum_cf3_tla_s_60khz = spectrum_cf3_tla_s_60khz / max(spectrum_cf3_tla_s_60khz(3500:4000));
% Rac-TFLA
spectrum_cf3_tla_rac_14khz = spectrum_cf3_tla_rac_14khz / max(spectrum_cf3_tla_rac_14khz(3500:4000));
spectrum_cf3_tla_rac_17p5khz = spectrum_cf3_tla_rac_17p5khz / max(spectrum_cf3_tla_rac_17p5khz(3500:4000));
spectrum_cf3_tla_rac_22khz = spectrum_cf3_tla_rac_22khz / max(spectrum_cf3_tla_rac_22khz(3500:4000));
spectrum_cf3_tla_rac_30khz = spectrum_cf3_tla_rac_30khz / max(spectrum_cf3_tla_rac_30khz(3500:4000));
spectrum_cf3_tla_rac_40khz = spectrum_cf3_tla_rac_40khz / max(spectrum_cf3_tla_rac_40khz(3500:4000));
spectrum_cf3_tla_rac_50khz = spectrum_cf3_tla_rac_50khz / max(spectrum_cf3_tla_rac_50khz(3500:4000));
spectrum_cf3_tla_rac_60khz = spectrum_cf3_tla_rac_60khz / max(spectrum_cf3_tla_rac_60khz(3500:4000));

% Make a nice plot of the spectra we have
spectra = {spectrum_cf3_tla_s_14khz, spectrum_cf3_tla_s_17p5khz, spectrum_cf3_tla_s_22khz, spectrum_cf3_tla_s_30khz, spectrum_cf3_tla_s_40khz, spectrum_cf3_tla_s_50khz, ...
    spectrum_cf3_tla_s_60khz, spectrum_cf3_tla_rac_14khz, spectrum_cf3_tla_rac_17p5khz, spectrum_cf3_tla_rac_22khz, spectrum_cf3_tla_rac_30khz, spectrum_cf3_tla_rac_40khz, ...
    spectrum_cf3_tla_rac_50khz,spectrum_cf3_tla_rac_60khz};


%% Variables section
% Bloch-McConnell variables
k_ex = 100; % extracted from chisquare fit in order: 286,286, 184, 184, 286, 184
T_2 = 0.01; % extracted from chisquare fit in order: 0.01 for all
J_cf = 280; % from experimental values (@Igor)

% NMR variables
t = linspace(0, 0.1, 8192); % time domain size (s) need to match experiment
DW = 5e-6; % dwell time
SW = 1/(2*DW); % spectral width (Hz)
sfrq = 700; % spectrometer freq (MHz)
offset = 9200; %offset (Hz) to align simulations and experiments ogni 100 scala di 2
NP = length(zeros(size(t))); % Length of the signal
f = linspace(-SW, SW, NP); % Frequency axis
f2 = linspace(-SW+offset, SW+offset, NP);
aqt = NP / (2*SW); % acquisition time of simulated spectrum

%% Simulate signals
% Modified Bloch equtions
L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
     k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit

L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
      3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
      0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
      0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];

% Compute and simplify the matrix exponential to get the time evolution
U1 = zeros(2, 2, length(t)); % Preallocate U1 matrix
for i = 1:length(t)
    U1(:,:,i) = expm(L1 * abs(t(i)));
end

U2 = zeros(4, 4, length(t)); % Preallocate U2 matrix
for i = 1:length(t)
    U2(:,:,i) = expm(L2 * abs(t(i)));
end

% Get the signals
vec_u1 = [2; 2];
vec_u2 = ones(4, 1);

% Compute signal 1
signal1 = zeros(size(t));
for i = 1:length(t)
    signal1(i) = sum(sum(U1(:,:,i) .* vec_u1));
end

% Compute signal2
signal2 = zeros(size(t));
for i = 1:length(t)
    signal2(i) = sum(sum(U2(:,:,i) .* vec_u2));
end

% Fourier transform and shift, see property of fft in matlab
spectrum1 = fft(signal1);
spectrum1 = fftshift(spectrum1);

spectrum2 = fft(signal2);
spectrum2 = fftshift(spectrum2);

total_spectrum = spectrum1 + spectrum2;

%% Spectral normalization section
real_part = real(total_spectrum);
imaginary_part = imag(total_spectrum);

% Normalize the real and imaginary parts separately
max_real = max(real_part(:));
min_real = min(real_part(:));
normalized_real = (real_part - min_real) / (max_real - min_real);

max_imaginary = max(imaginary_part(:));
min_imaginary = min(imaginary_part(:));
normalized_imaginary = (imaginary_part - min_imaginary) / (max_imaginary - min_imaginary);

% Combine the normalized real and imaginary parts to get the normalized spectrum
normalized_spectrum = normalized_real; %+ 1i * normalized_imaginary;
normalized_spectrum = normalized_spectrum';

% Interpolate the shifted spectrum
shifted_spectrum = interp1(f, normalized_spectrum, f2, 'linear', 0);

%% Produce plots
% plot normalized spectra on each other
figure(1); clf; hold on; 
plot(shifted_spectrum)
plot(normalized_spectrum)
plot(spectrum_cf3_tla_rac_60khz)
xlabel('Hz')
legend('Simulated - offset', 'Simulated - original', 'Experimental')












