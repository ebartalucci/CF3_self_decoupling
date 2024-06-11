% Chisquare fit to extract K_ex and T_2
% Author: Ettore Bartalucci, RWTH Aachen
% Scripts for Bloch-McConnel from Matthias Ernst, ETH Zurich
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 19.04.24
% Project: CF3 self decoupling

clear all;
tic;

si=8192;
name= 'spectra/TLA_S_14khz_exp_200_161023/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_14khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/TLA_S_30khz_exp_14_121023/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_30khz = fread(fid,si,'int');
fclose(fid);

si=131072;
name= 'spectra/TLA_S_60khz_exp_15_121023/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_s_60khz = fread(fid,si,'int');
fclose(fid);

si=262144;
name= 'spectra/TLA_rac_14khz_exp_10_221123/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_14khz = fread(fid,si,'int');
fclose(fid);

si=8192;
name= 'spectra/TLA_rac_30khz_exp_10_211123/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_30khz = fread(fid,si,'int');
fclose(fid);

si=131072;
name= 'spectra/TLA_rac_60khz_exp_13_211123/pdata/1/1r';
fid = fopen(name,'r','l');
spectrum_cf3_tla_rac_60khz = fread(fid,si,'int');
fclose(fid);

% Normalize experimental spectral intensity of CF3 peak between 0 and 1
spectrum_cf3_tla_s_14khz = spectrum_cf3_tla_s_14khz / max(spectrum_cf3_tla_s_14khz(3294:3794));
spectrum_cf3_tla_s_30khz = spectrum_cf3_tla_s_30khz / max(spectrum_cf3_tla_s_30khz(3500:4000));
spectrum_cf3_tla_s_60khz = spectrum_cf3_tla_s_60khz / max(spectrum_cf3_tla_s_60khz(55000:65000));
spectrum_cf3_tla_rac_14khz = spectrum_cf3_tla_rac_14khz / max(spectrum_cf3_tla_rac_14khz(114922:124922));
spectrum_cf3_tla_rac_30khz = spectrum_cf3_tla_rac_30khz / max(spectrum_cf3_tla_rac_30khz(3470:3970));
spectrum_cf3_tla_rac_60khz = spectrum_cf3_tla_rac_60khz / max(spectrum_cf3_tla_rac_60khz(54535:64535));

% Make a nice plot of the spectra we have
spectra = {spectrum_cf3_tla_s_14khz, spectrum_cf3_tla_s_30khz, spectrum_cf3_tla_s_60khz, ...
           spectrum_cf3_tla_rac_14khz, spectrum_cf3_tla_rac_30khz, spectrum_cf3_tla_rac_60khz};

subplot_titles = {'S-TFLA 14 kHz','S-TFLA 30 kHz','S-TFLA 60 kHz',...
                  'rac-TFLA 14 kHz','rac-TFLA 30 kHz','rac-TFLA 60 kHz'};

% figure(999);
% for i = 1:6
%     subplot(3, 2, i);
%     plot(spectra{i});
%     title(subplot_titles{i});
%     ylim([0,1.2]);
% end

%% Simulations section
% Define range of k_ex and T_2
k_ex_values = linspace(1, 1000, 50); % change this to 50 or 100, takes long
T_2_values = linspace(0.01, 0.1, 50); % change this to 50 or 100, takes long
J_cf = 280; % exp value

% Define different t values that fits the experimental sizes
t_values = {linspace(0, 8191*1e-5, 8192), linspace(0, 8191*1e-5, 8192), ...
            linspace(0, 131071*1e-5, 131072), linspace(0, 262143*1e-5, 262144), ...
            linspace(0, 8191*1e-5, 8192), linspace(0, 131071*1e-5, 131072)};

% Printing file
fileID = fopen('chisquare_mins.txt','w');

% Loop over each TLA spectrum
for t_index = 1:length(t_values)

    % Select the current spectrum based on the t_index, dont change order
    switch t_index
        case 1
            current_spectrum = spectrum_cf3_tla_s_14khz;
            offset = 6750; %offset (Hz) to align simulations and experiments
        case 2
            current_spectrum = spectrum_cf3_tla_s_30khz;
            offset = 4250; %offset (Hz) to align simulations and experiments
        case 3
            current_spectrum = spectrum_cf3_tla_s_60khz;
            offset = 4250;
        case 4
            current_spectrum = spectrum_cf3_tla_rac_14khz;
            offset = 4250;
        case 5
            current_spectrum = spectrum_cf3_tla_rac_30khz;
            offset = 4600; %offset (Hz) to align simulations and experiments 
        case 6
            current_spectrum = spectrum_cf3_tla_rac_60khz;
            offset = 4600;
    end
    
    % NMR variables
    t = t_values{t_index}; % time domain size (s) need to match experiment
    DW = t(2) - t(1); % dwell time
    SW = 1/(DW); % spectral width (Hz)
    NP = length(zeros(size(t))); % Length of the signal
      
    % Store chi-square statistics
    chi_square_matrix = zeros(length(k_ex_values), length(T_2_values));

    % Compute the  frequency axis
    f = linspace(-SW/2, SW/2, NP); % Frequency axis
    f2 = linspace(-SW/2+offset, SW/2+offset, NP); % Shifted axis by offset

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
    
    %%%%%%%%%%%%%%%%%%%% END MATTHIAS ERNST, ETH ZURICH %%%%%%%%%%%%%%%%%%%
            
            final_spectrum = real(final_spectrum(:))/max(real(final_spectrum));
%            plot(f,final_spectrum,f2,current_spectrum)
            [a1,b1a] = min(abs(f+3000));
            [a1,b1b] = min(abs(f-3000));
            [a1,b2a] = min(abs(f2+3000));
            b2b=b2a+b1b-b1a;
%    		plot(f(b1a:b1b),final_spectrum(b1a:b1b),f2(b2a:b2b),current_spectrum(b2a:b2b))
            chi_square_matrix(k, T) = sum((final_spectrum(b1a:b1b)-current_spectrum(b2a:b2b)).^2);
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
    disp(['Minimal Chi-square for TLA Spectrum ', num2str(t_index), ':']);
    disp(['T2 Value: ', num2str(min_T_2)]);
    disp(['k_ex Value: ', num2str(min_k_ex)]);
    disp(['Chi-square Value: ', num2str(min_chi_square)]);
    disp(' ');
    
    % Write the results to the text file
    fprintf(fileID, 'Minimal Chi-square for TLA Spectrum %d:\n', t_index);
    fprintf(fileID, 'T2 Value: %f\n', min_T_2);
    fprintf(fileID, 'k_ex Value: %f\n', min_k_ex);
    fprintf(fileID, 'Chi-square Value: %f\n\n', min_chi_square);

    % Normalize the chi-square matrix to the range [0, 1]
    normalized_chi_square_matrix = chi_square_matrix / max(max(chi_square_matrix));

    % Plot the chi-square statistics as a heatmap with a specified colormap
    subplot(2, 3, t_index); 
    imagesc(T_2_values, k_ex_values, normalized_chi_square_matrix);
    colormap(jet); 
    colorbar; 
    xlabel('T2');
    ylabel('k_ex');
    title('Chi-square'); % Need to set title using spectrum_titles array
end

% Close the text file
fclose(fileID);

runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%