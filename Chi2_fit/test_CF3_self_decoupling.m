% Simulations of Self-decoupling as exchange process
% Originally written by Matthias Ernst, ETH Zurich
% Matlab adaptation: Ettore Bartalucci, RWTH Aachen
% First draft: Aachen, 07.03.24
% Last update: Aachen, 04.04.24
% Project: CF3 self decoupling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TEST: load experimental data                            %%%%%%%%%

% Open spectrum using RBNMR function and plot to see if everything alright
NMR_data = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');

% Select region in the CF3 region of the spectrum (115ppm-135ppm)
w = NMR_data.XAxis;
spectrum = NMR_data.Data;

w_cf3 = w(3294:3794);
spectrum_cf3 = spectrum(3294:3794);
spectrum_cf3 = spectrum_cf3 / max(spectrum_cf3);

plot(spectrum_cf3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TEST: Simulate one spectrum for k_ex=0.001 and T_2=0.03 %%%%%%%%%

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

figure(2);
plot(f, normalized_spectrum)
runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


