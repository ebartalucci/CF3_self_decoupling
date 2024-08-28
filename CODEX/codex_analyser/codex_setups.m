%% AVOID RELOADING
SkipVars = {'FOLDER'};

for i = 1:length(SkipVars)
    varName = SkipVars{i};
    if evalin('base', ['exist(''', varName, ''', ''var'')'])
        disp([varName ' already exists in the workspace, skipping...']);
    else
        disp([varName ' is not in the workspace, loading...']);
        FOLDER = rbnmr('D:\PhD\Data\NMR\ETBA\Bolm\CF3_coupling\190824_3p2mm_HX_CODEX_setup_TWC_EB_NC');
    end
end

Spectra = FOLDER(17:36);

SB1 = 22000:26000;
SB2 = 38000:41000;
NB = 9000:11000;

BinS0 = 1; % 1 if S0 is the second experiment of the pair, 0 if its the first

Offset = 1.8e11; % for stacked plots

Integ_1_S = []; 
Integ_1_S0 = []; 
Integ_2_S = []; 
Integ_2_S0 = []; 

%% SPLIT INTO S AND S0
if BinS0 == 1
    Spectrum_S = Spectra(1:2:end); % odd indexes
    Spectrum_S0 = Spectra(2:2:end); % even indexes
elseif BinS0 == 0
    Spectrum_S = Spectra(2:2:end); % even indexes
    Spectrum_S0 = Spectra(1:2:end); % odd indexes
end

%% GET SPECTRA
for i = 1:length(Spectrum_S)
    StructArrayS = Spectrum_S{i};
    StructArrayS0 = Spectrum_S0{i};
    for j = 1:length(StructArrayS)
        S_spectrum = StructArrayS(j).Data;
        S0_spectrum = StructArrayS0(j).Data;
        
        I_S_1 = sum(S_spectrum(SB1), 2); % integrals
        I_S_1 = sum(I_S_1(:));

        I_S0_1 = sum(S0_spectrum(SB1), 2); % integrals
        I_S0_1 = sum(I_S0_1(:));
        
        I_S_2 = sum(S_spectrum(SB2), 2); % integrals
        I_S_2 = sum(I_S_2(:));

        I_S0_2 = sum(S0_spectrum(SB2), 2); % integrals
        I_S0_2 = sum(I_S0_2(:));

        % Stack the spectra by adding an offset
        S_spectrum_1 = S_spectrum(SB1) + (i-1)*Offset;
        S0_spectrum_1 = S0_spectrum(SB1) + (i-1)*Offset;

        S_spectrum_2 = S_spectrum(SB2) + (i-1)*Offset;
        S0_spectrum_2 = S0_spectrum(SB2) + (i-1)*Offset;
        
        Integ_1_S(end+1) = I_S_1; 
        Integ_1_S0(end+1) = I_S0_1;

        Integ_2_S(end+1) = I_S_2; 
        Integ_2_S0(end+1) = I_S0_2;

        % plot spectra
        figure(2)
        subplot(2,2,1)
        plot(S_spectrum_1)
        title('Glycine peak 1 1H-13C CP S CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([1200 4000])
        ylim([-1e11 2.3e12])
        axis square
        
        hold on;

        subplot(2,2,2)
        plot(S0_spectrum_1)
        title('Glycine peak 1 1H-13C CP S0 CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([1200 4000])
        ylim([-1e11 2.3e12])
        axis square

        hold on;

        subplot(2,2,3)
        plot(S_spectrum_2)
        title('Glycine peak 2 1H-13C CP S CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e11 2.5e12])
        axis square
        
        hold on;

        subplot(2,2,4)
        plot(S0_spectrum_2)
        title('Glycine peak 2 1H-13C CP S0 CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e11 2.5e12])
        axis square

        hold on;

    end
end

set(gcf, 'PaperOrientation', 'landscape'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 29.7 21]); 
set(gcf, 'PaperSize', [29.7 21]); % A4 size in centimeters (landscape)
print(gcf, 'Setup_CODEX_plots', '-dpdf', '-r300');

%% NORMALIZE S/S0 AND PLOT
n = length(Integ_1_S);
if mod(n, 2) ~= 0
    Integ_1_S = Integ_1_S(1:end-1); % Remove the last element if odd
    Integ_1_S0 = Integ_1_S0(1:end-1); % Remove the last element if odd
end

q = length(Integ_2_S);
if mod(q, 2) ~= 0
    Integ_2_S = Integ_2_S(1:end-1); % Remove the last element if odd
    Integ_2_S0 = Integ_2_S0(1:end-1); % Remove the last element if odd
end

norm_CODEX_1 = Integ_1_S ./ Integ_1_S0;
norm_CODEX_2 = Integ_2_S ./ Integ_2_S0;

% Plot the results
figure(3);
subplot(2,1,1)
plot(norm_CODEX_1, '-o');
title('rac-TFLA CODEX Evolution')
xlabel('Mixing Time')
ylabel('S / S0')
grid on;
subplot(2,1,2)
plot(norm_CODEX_2, '-o');
title('S-TFLA CODEX Evolution')
xlabel('Mixing Time')
ylabel('S / S0')
grid on;

set(gcf, 'PaperOrientation', 'landscape'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 29.7 21]); 
set(gcf, 'PaperSize', [29.7 21]); % A4 size in centimeters (landscape)
print(gcf, 'Setup_CODEX_decay_plots', '-dpdf', '-r300');
