% Chisquare fit to extract K_ex and T_2
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 08.04.24
% Project: CF3 self decoupling

tic;

%% Experimental section
% Optical pure TLA @ 14kHz, 30kHz, and 60kHz
tla_s_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');
tla_s_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');
tla_s_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_60khz_exp_15_121023\pdata\1');
% Racemic @ 14kHz, 30kHz, and 60kHz
tla_rac_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_14khz_exp_10_221123\pdata\1');
tla_rac_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_30khz_exp_10_211123\pdata\1');
tla_rac_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_60khz_exp_13_211123\pdata\1');

% Select CF3 region of the spectrum (115ppm-135ppm), centered on max signal
% Extract CF3 spectra for each TLA spectrum
spectrum_cf3_tla_s_14khz = tla_s_14khz.Data(3294:3794); % here are 500 points from middle of spectrum
spectrum_cf3_tla_s_30khz = tla_s_30khz.Data(3500:4000);
% spectrum_cf3_tla_s_60khz = tla_s_60khz.Data(4750:7250);
% spectrum_cf3_tla_rac_14khz = tla_rac_14khz.Data(118672:121172);
spectrum_cf3_tla_rac_30khz = tla_rac_30khz.Data(3470:3970);
% spectrum_cf3_tla_rac_60khz = tla_rac_60khz.Data(58285:60785);
% 
% figure(1212);hold on;
% plot(spectrum_cf3_tla_s_14khz)
% plot(spectrum_cf3_tla_s_30khz)
% plot(spectrum_cf3_tla_rac_30khz)
% legend();
% hold off;


%% Simulations section
% Define range of k_ex and T_2
k_ex_values = linspace(1, 1000, 10); % change this to 50
T_2_values = linspace(0.01, 0.1, 10); % change this to 50
J_cf = 280; 
t = linspace(0, 0.1, 501);
Fs = 1 / (t(2) - t(1)); % Sampling frequency (Hz)

% printing file
fileID = fopen('chisquare_mins.txt','w');

% Loop over each TLA spectrum
for spectrum_index = 1:3
    
    % store chi-square statistics
    chi_square_matrix = zeros(length(k_ex_values), length(T_2_values));

    % Select the current TLA spectrum
    switch spectrum_index
        case 1
            current_spectrum = spectrum_cf3_tla_s_14khz;
        case 2
            current_spectrum = spectrum_cf3_tla_s_30khz;
        % case 3
        %     current_spectrum = spectrum_cf3_tla_s_60khz;
        % case 4
        %     current_spectrum = spectrum_cf3_tla_rac_14khz;
        case 3
            current_spectrum = spectrum_cf3_tla_rac_30khz;
        % case 6
        %     current_spectrum = spectrum_cf3_tla_rac_60khz;
    end

    %%%%%%%%%%%%%%%%%%% FROM MATTHIAS ERNST, ETH ZURICH %%%%%%%%%%%%%%%%%%%
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
    
            %%%%%%%%%%%%%%%%%%%%% END MATTHIAS ERNST %%%%%%%%%%%%%%%%%%%%%%
            
            % Normalize final spectrum by splitting real and imaginary and
            % recombining them at the end
            real_final_spectrum = real(final_spectrum);
            imaginary_final_spectrum = imag(final_spectrum);
            
            % Normalize the real and imaginary parts separately
            max_real = max(real_final_spectrum(:));
            min_real = min(real_final_spectrum(:));
            normalized_real_final_spectrum = (real_final_spectrum - min_real) / (max_real - min_real);
            
            max_imaginary = max(imaginary_final_spectrum(:));
            min_imaginary = min(imaginary_final_spectrum(:));
            normalized_imaginary_final_spectrum = (imaginary_final_spectrum - min_imaginary) / (max_imaginary - min_imaginary);
            
            % Combine the normalized real and imaginary parts to get the normalized spectrum
            norm_final_spectrum = normalized_real_final_spectrum + 1i * normalized_imaginary_final_spectrum;      

            % Compute the  frequency axis
            N = length(signal1); % Length of the signal
            f = linspace(-Fs/2, Fs/2, N); % Frequency axis

            % normalize experimental spectrum
            norm_current_spectrum = current_spectrum / max(current_spectrum);
            norm_current_spectrum = norm_current_spectrum';

            % Compute chi-square statistic between experimental and simulated spectra
            chi_square = sum((norm_current_spectrum - norm_final_spectrum).^2 ./ norm_final_spectrum);

            % Store chi-square statistic in the matrix
            chi_square_matrix(k, T) = real(chi_square);
        end
    end

    %% Plot chi-square statistics

    % Compute the indices of the minimum chi-square value
    [min_chi_square, min_index] = min(chi_square_matrix(:));
    [min_k_ex_index, min_T_2_index] = ind2sub(size(chi_square_matrix), min_index);
    
    % Retrieve the corresponding k_ex and T_2 values
    min_k_ex = k_ex_values(min_k_ex_index);
    min_T_2 = T_2_values(min_T_2_index);
    
    % Print the results to the console
    disp(['Minimal Chi-square for TLA Spectrum ', num2str(spectrum_index), ':']);
    disp(['T2 Value: ', num2str(min_T_2)]);
    disp(['k_ex Value: ', num2str(min_k_ex)]);
    disp(['Chi-square Value: ', num2str(min_chi_square)]);
    disp(' ');
    
    % Write the results to the text file
    fprintf(fileID, 'Minimal Chi-square for TLA Spectrum %d:\n', spectrum_index);
    fprintf(fileID, 'T2 Value: %f\n', min_T_2);
    fprintf(fileID, 'k_ex Value: %f\n', min_k_ex);
    fprintf(fileID, 'Chi-square Value: %f\n\n', min_chi_square);

    % Normalize the chi-square matrix to the range [0, 1]
    normalized_chi_square_matrix = (chi_square_matrix - min(chi_square_matrix(:))) / (max(chi_square_matrix(:)) - min(chi_square_matrix(:)));

    % Plot the chi-square statistics as a heatmap with a specified colormap
    subplot(2, 3, spectrum_index); 
    imagesc(T_2_values, k_ex_values, normalized_chi_square_matrix);
    colormap(jet); 
    colorbar; 
    xlabel('T2');
    ylabel('k_ex');
    title('Chi-square'); % Set title using spectrum_titles array
end

% Close the text file
fclose(fileID);

runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%