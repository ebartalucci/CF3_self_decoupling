%% Experimental section
clear all;

% Load the experimental spectra
% Optical pure TLA @ 14kHz, 30kHz, and 60kHz
tla_s_14khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');
tla_s_30khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');
tla_s_60khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_60khz_exp_15_121023\pdata\1');
% Racemic @ 14kHz, 30kHz, and 60kHz
tla_rac_14khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_14khz_exp_10_221123\pdata\1');
tla_rac_30khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_30khz_exp_10_211123\pdata\1');
tla_rac_60khz = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_60khz_exp_13_211123\pdata\1');

% Select CF3 region of the spectrum (115ppm-135ppm)
% For each spectrum, extract the CF3 region
% You can define this based on your specific data
cf3_region_indices = 3294:3794; % Indices corresponding to CF3 region

% Extract CF3 spectra for each TLA spectrum
spectrum_cf3_tla_s_14khz = tla_s_14khz.Data(3294:3794); % here are 250 points from middle of spectrum
spectrum_cf3_tla_s_30khz = tla_s_30khz.Data(3500:4000); % here are 250 points from middle of spectrum
spectrum_cf3_tla_s_60khz = tla_s_60khz.Data(55000:65000); % here are 5000 points from middle of spectrum
spectrum_cf3_tla_rac_14khz = tla_rac_14khz.Data(114922:124922); % here are 5000 points from middle of spectrum
spectrum_cf3_tla_rac_30khz = tla_rac_30khz.Data(3470:3970); % here are 250 points from middle of spectrum
spectrum_cf3_tla_rac_60khz = tla_rac_60khz.Data(54535:64535); % here are 5000 points from middle of spectrum


spectrum_cf3_tla_rac_60khz = spectrum_cf3_tla_rac_60khz / max(spectrum_cf3_tla_rac_60khz);


% Initialize variables
k_ex = 150;
T_2 = 0.1;
J_cf = 280; 
t = linspace(0, 0.1, 501);
Fs = 1/(t(2) - t(1)); % Sampling frequency (Hz)


% Simulate signals
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

% Fourier transform and shift 
spectrum1 = fft(signal1);
spectrum1 = fftshift(spectrum1);

spectrum2 = fft(signal2);
spectrum2 = fftshift(spectrum2);

% Compute the  frequency axis
N = length(signal1); % Length of the signal
f = linspace(-Fs/2, Fs/2, N); % Frequency axis

total_spectrum = spectrum1 + spectrum2;


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
normalized_spectrum = normalized_real + 1i * normalized_imaginary;
normalized_spectrum = normalized_spectrum';
spectrum_cf3_tla_s_30khz = spectrum_cf3_tla_s_30khz / max(spectrum_cf3_tla_s_30khz);

figure(999);hold on;
plot(f, normalized_spectrum)
plot(f, spectrum_cf3_tla_s_30khz)

% Calculate chi-square without division by zero
real_spectrum_cf3_tla_s_30khz = real(spectrum_cf3_tla_s_30khz);
real_normalized_spectrum = real(normalized_spectrum);
non_zero_indices = real_normalized_spectrum ~= 0;

chisquare = sum((real_spectrum_cf3_tla_s_30khz(non_zero_indices) - real_normalized_spectrum(non_zero_indices)).^2 ./ real_normalized_spectrum(non_zero_indices));






runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);

