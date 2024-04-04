% Chisquare fit to extract K_ex and T_2
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 04.04.24
% Project: CF3 self decoupling

%% Simulations section
% For the simulations: k_ex=(1, 1000), T_2=(1, 5) 
tic;

% Define range of k_ex and T_2
k_ex_values = linspace(1, 1000, 10); 
T_2_values = linspace(0.01, 0.1, 10); 

% Plotting for fixed T_2 = 3
for k = 1:length(k_ex_values)
    % variables
    k_ex = k_ex_values(k);
    T_2 = 3/100;
    J_cf = 280; 
    t = linspace(0, 0.1, 1000);
    Fs = 1/(t(2) - t(1)); % Sampling frequency (Hz)

    % Modified Bloch equtions
    L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
         k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit
    
    L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
          3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
          0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
          0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];

    % Compute the matrix exponentials to get the time evolutions
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

    final_spectrum = spectrum1+spectrum2;
    
    % Compute the  frequency axis
    N = length(signal1); % Length of the signal
    f = linspace(-Fs/2, Fs/2, N); % Frequency axis
   
    % figure;
    % plot(f, abs(final_spectrum), 'LineWidth', 1.5);
    % xlabel('Frequency (Hz)', 'FontSize', 12);
    % title(['Final Spectrum for(k_{ex} = ', num2str(k_ex), ', T_2 = 3)'], 'FontSize', 12);
end

% Plotting for fixed k_ex = 100
for T = 1:length(T_2_values)
    T_2 = T_2_values(T);
    k_ex = 100;
    J_cf = 280; 
    t = linspace(0, 0.1, 1000);
    Fs = 1/(t(2) - t(1)); % Sampling frequency (Hz)

    % Modified Bloch equtions
    L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
         k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit
    
    L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
          3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
          0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
          0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];

    % Compute the matrix exponentials to get the time evolutions
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

    final_spectrum = spectrum1+spectrum2;
    
    % Compute the  frequency axis
    N = length(signal1); % Length of the signal
    f = linspace(-Fs/2, Fs/2, N); % Frequency axis
   
    % figure;
    % plot(f, abs(final_spectrum), 'LineWidth', 1.5);
    % xlabel('Frequency (Hz)', 'FontSize', 12);
    % title(['Plot of Signal1 (k_{ex} = 100, T_2 = ', num2str(T_2), ')'], 'FontSize', 12);

end

% Initialize final spectrum matrix
final_spectrum_matrix = zeros(length(k_ex_values), length(T_2_values));

% Compute final spectrum for all combinations of k_ex and T_2 values
for k = 1:length(k_ex_values)
    for T = 1:length(T_2_values)
        k_ex = k_ex_values(k);
        T_2 = T_2_values(T);

        % Modified Bloch equtions
        L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
             k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit
        
        L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
              3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
              0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
              0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];
    
        % Compute the matrix exponentials to get the time evolutions
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
    
        % Compute final spectrum
        final_spectrum = spectrum1 + spectrum2;
        final_spectrum_matrix(k, T) = max(abs(final_spectrum)); % Store max amplitude
    end
end

% Plot 2D contour map
[X, Y] = meshgrid(k_ex_values, T_2_values);

figure;
contourf(X, Y, final_spectrum_matrix', 'LineWidth', 1.5);
colorbar;
xlabel('k_{ex}', 'FontSize', 12);
ylabel('T_2', 'FontSize', 12);
title('2D Contour Map of Final Spectrum', 'FontSize', 12);

runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%