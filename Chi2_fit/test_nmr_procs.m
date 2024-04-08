%% Experimental section
% Load the experimental spectra
% Optical pure @ 14kHz, 30kHz, and 60kHz
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
spectrum_cf3_tla_s_14khz = tla_s_14khz.Data(cf3_region_indices);
spectrum_cf3_tla_s_30khz = tla_s_30khz.Data(cf3_region_indices);
spectrum_cf3_tla_s_60khz = tla_s_60khz.Data(cf3_region_indices);
spectrum_cf3_tla_rac_14khz = tla_rac_14khz.Data(cf3_region_indices);
spectrum_cf3_tla_rac_30khz = tla_rac_30khz.Data(cf3_region_indices);
spectrum_cf3_tla_rac_60khz = tla_rac_60khz.Data(cf3_region_indices);

%% Simulations section
% Define range of k_ex and T_2
k_ex_values = linspace(1, 1000, 10); 
T_2_values = linspace(0.01, 0.1, 10);
J_cf = 280; 
t = linspace(0, 0.1, 501);
Fs = 1 / (t(2) - t(1)); % Sampling frequency (Hz)

% Loop over each TLA spectrum
for spectrum_index = 1:6
    % Initialize matrix to store chi-square statistics
    chi_square_matrix = zeros(length(k_ex_values), length(T_2_values));

    % Select the current TLA spectrum
    switch spectrum_index
        case 1
            current_spectrum = spectrum_cf3_tla_s_14khz;
        case 2
            current_spectrum = spectrum_cf3_tla_s_30khz;
        case 3
            current_spectrum = spectrum_cf3_tla_s_60khz;
        case 4
            current_spectrum = spectrum_cf3_tla_rac_14khz;
        case 5
            current_spectrum = spectrum_cf3_tla_rac_30khz;
        case 6
            current_spectrum = spectrum_cf3_tla_rac_60khz;
    end

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
    
            % Compute the  frequency axis
            N = length(signal1); % Length of the signal
            f = linspace(-Fs/2, Fs/2, N); % Frequency axis

            % Compute chi-square statistic between experimental and simulated spectra
            % Use the corresponding CF3 spectrum for the current TLA spectrum
            chi_square = sum((current_spectrum' - final_spectrum).^2 ./ final_spectrum);

            % Store chi-square statistic in the matrix
            chi_square_matrix(k, T) = real(chi_square);
        end
    end

    %% Plot chi-square statistics
    % Normalize the chi-square matrix to the range [0, 1]
    normalized_chi_square_matrix = (chi_square_matrix - min(chi_square_matrix(:))) / (max(chi_square_matrix(:)) - min(chi_square_matrix(:)));

    % Plot the chi-square statistics as a heatmap with a specified colormap
    subplot(2, 3, spectrum_index); % Adjust subplot position based on your preference
    imagesc(T_2_values, k_ex_values, normalized_chi_square_matrix);
    colormap(jet); % Specify the colormap
    colorbar; % Display colorbar
    xlabel('T2 Values');
    ylabel('k_ex Values');
    title(sprintf('Chi-square Statistics Heatmap - TLA Spectrum %d', spectrum_index));
end
