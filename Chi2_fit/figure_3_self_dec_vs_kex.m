% Generate figure for CF3 self-decoupling at various exchange rates
% Simulations of Self-decoupling as exchange process
% Originally written by Matthias Ernst, ETH Zurich
% Matlab adaptation: Ettore Bartalucci, RWTH Aachen
% Project: CF3 self decoupling

%% Variables section
% Bloch-McConnell variables
k_ex_values = [1, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500]; 
T_2 = 3/100; 
lw = 1/(pi * T_2);
disp(lw)
J_cf = 280; 

% NMR variables
t = linspace(0, 0.1, 200); % time domain size (s) need to match experiment
Fs = 1/(t(2) - t(1)); % Frequency axis
% Compute the  frequency axis
N = length(t); % Length of the signal
f = linspace(-Fs/2, Fs/2, N); % Frequency axis

%% Simulate signals
for k = 1:length(k_ex_values)

    k_ex = k_ex_values(k)

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
    
    %% Produce plots
    % plot normalized spectra on each other
    subplot(3, 4, k)
    plot(f, normalized_spectrum)
    ylim([-0.05, 1.1])
    switch k
        case 1
            title('k_{ex} = 1 s^{-1}')
        case 2
            title('k_{ex} = 100 s^{-1}')
        case 3
            title('k_{ex} = 200 s^{-1}')
        case 4
            title('k_{ex} = 300 s^{-1}')
        case 5
            title('k_{ex} = 500 s^{-1}')
        case 6
            title('k_{ex} = 750 s^{-1}')
        case 7
            title('k_{ex} = 1000 s^{-1}')
        case 8
            title('k_{ex} = 1500 s^{-1}')
        case 9
            title('k_{ex} = 2000 s^{-1}')
        case 10
            title('k_{ex} = 3000 s^{-1}')
        case 11
            title('k_{ex} = 5000 s^{-1}')
        case 12
            title('k_{ex} = 7500 s^{-1}')
    end
    xlabel('\nu [kHz]')
end