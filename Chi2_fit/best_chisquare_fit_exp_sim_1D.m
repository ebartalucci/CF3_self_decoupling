% Compare experiments and simulations with best T2 and k_ex
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 18.04.24
% Project: CF3 self decoupling

%% Experimental section
clear all;

% Load the experimental spectra
% Optical pure TLA @ 14kHz, 30kHz, and 60kHz
tla_s_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');
tla_s_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');
tla_s_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_60khz_exp_15_121023\pdata\1');
% Racemic @ 14kHz, 30kHz, and 60kHz
tla_rac_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_14khz_exp_10_221123\pdata\1');
tla_rac_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_30khz_exp_10_211123\pdata\1');
tla_rac_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_60khz_exp_13_211123\pdata\1');

% tla_s_14khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');
% tla_s_30khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');
% tla_s_60khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_60khz_exp_15_121023\pdata\1');
% % Racemic @ 14kHz, 30kHz, and 60kHz
% tla_rac_14khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_14khz_exp_10_221123\pdata\1');
% tla_rac_30khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_30khz_exp_10_211123\pdata\1');
% tla_rac_60khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_60khz_exp_13_211123\pdata\1');


% Extract CF3 spectra for each TLA spectrum
spectrum_cf3_tla_s_14khz = tla_s_14khz.Data(3294:3794); % here are 250 points from middle of spectrum
spectrum_cf3_tla_s_30khz = tla_s_30khz.Data(3500:4000); % here are 250 points from middle of spectrum
spectrum_cf3_tla_s_60khz = tla_s_60khz.Data(55000:65000); % here are 5000 points from middle of spectrum
spectrum_cf3_tla_rac_14khz = tla_rac_14khz.Data(114922:124922); % here are 5000 points from middle of spectrum
spectrum_cf3_tla_rac_30khz = tla_rac_30khz.Data;%(3470:3970); % here are 250 points from middle of spectrum
spectrum_cf3_tla_rac_60khz = tla_rac_60khz.Data(54535:64535); % here are 5000 points from middle of spectrum

% Normalize experimental spectral intensity between 0 and 1
spectrum_cf3_tla_s_14khz = spectrum_cf3_tla_s_14khz / max(spectrum_cf3_tla_s_14khz);
spectrum_cf3_tla_s_30khz = spectrum_cf3_tla_s_30khz / max(spectrum_cf3_tla_s_30khz);
spectrum_cf3_tla_s_60khz = spectrum_cf3_tla_s_60khz / max(spectrum_cf3_tla_s_60khz);
spectrum_cf3_tla_rac_14khz = spectrum_cf3_tla_rac_14khz / max(spectrum_cf3_tla_rac_14khz);
spectrum_cf3_tla_rac_30khz = spectrum_cf3_tla_rac_30khz / max(spectrum_cf3_tla_rac_30khz(3470:3970));
spectrum_cf3_tla_rac_60khz = spectrum_cf3_tla_rac_60khz / max(spectrum_cf3_tla_rac_60khz);

%% Variables section
% Bloch-McConnell variables
k_ex = 60; % extracted from chisquare fit in order: 286,286, 184, 184, 286, 184
T_2 = 0.017; % extracted from chisquare fit in order: 0.01 for all
J_cf = 280; % from experimental values (@Igor)

% NMR variables
t = linspace(0, 0.1, 8192); % time domain size (s) need to match experiment
DW = 0.0000019; % dwell time
SW = 1/(2*DW); % spectral width (Hz)
sfrq = 700; % spectrometer freq (MHz)
offset = 24200; %offset (Hz) to align simulations and experiments ogni 100 scala di 2
NP = length(zeros(size(t))); % Length of the signal
f = linspace(-SW, SW, NP); % Frequency axis
f2 = linspace(-SW+offset, SW+offset, NP);
aqt = NP / (2*SW); % acquisition time of simulated spectrum
dig_res_sim_rac_30khz = 1 / aqt; %digital resolution, if everything right matches the exp one

%% Check spectral parameters
% Check digital resolution experimental spectra
dig_res_exp_rac_14khz = tla_rac_14khz.Acqus.SW_h / (tla_rac_14khz.Acqus.TD / 2);
dig_res_exp_rac_30khz = tla_rac_30khz.Acqus.SW_h / (tla_rac_30khz.Acqus.TD / 2);
dig_res_exp_rac_60khz = tla_rac_60khz.Acqus.SW_h / (tla_rac_60khz.Acqus.TD / 2);
dig_res_exp_s_14khz = tla_s_14khz.Acqus.SW_h / (tla_s_14khz.Acqus.TD / 2);
dig_res_exp_s_30khz = tla_s_30khz.Acqus.SW_h / (tla_s_30khz.Acqus.TD / 2);
dig_res_exp_s_60khz = tla_s_60khz.Acqus.SW_h / (tla_s_60khz.Acqus.TD / 2);
disp('Digital resolution of experimental spectra are all the same!')

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
plot(spectrum_cf3_tla_rac_30khz)
xlabel('Hz')
legend('Simulated - offset', 'Simulated - original', 'Experimental')
title('TLA-rac 30kHz')
























