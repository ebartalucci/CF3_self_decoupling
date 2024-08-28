%% AVOID RELOADING
SkipVars = {'FOLDER'};

for i = 1:length(SkipVars)
    varName = SkipVars{i};
    if evalin('base', ['exist(''', varName, ''', ''var'')'])
        disp([varName ' already exists in the workspace, skipping...']);
    else
        disp([varName ' is not in the workspace, loading...']);
        FOLDER = rbnmr('D:\PhD\Data\NMR\ETBA\Bolm\CF3_coupling\130824_3p2mm_19F_TLA_Static_and_CODEX');
    end
end

%% LOAD DATA
rac_TFLA = FOLDER(13:44);
S_TFLA = FOLDER(46:77);

SB = 7000:9400; % signal boundary
NB = 9500:10500; % noise boundary

BinS0 = 1; % 1 if S0 is the second experiment of the pair, 0 if its the first

Offset = 15000; % for stacked plots

Integ_rac_S = []; 
Integ_rac_S0 = []; 
Integ_S_S = []; 
Integ_S_S0 = []; 

%% SPLIT INTO S AND S0
if BinS0 == 1
    rac_TFLA_S = rac_TFLA(1:2:end); % odd indexes
    rac_TFLA_S0 = rac_TFLA(2:2:end); % even indexes
    S_TFLA_S = S_TFLA(1:2:end); % odd indexes
    S_TFLA_S0 = S_TFLA(2:2:end); % even indexes
elseif BinS0 == 0
    rac_TFLA_S = rac_TFLA(2:2:end); % even indexes
    rac_TFLA_S0 = rac_TFLA(1:2:end); % odd indexes
    S_TFLA_S = S_TFLA(2:2:end); % even indexes
    S_TFLA_S0 = S_TFLA(1:2:end); % odd indexes
end

%% GET SPECTRA
for i = 1:length(rac_TFLA_S)
    StructArrayS = rac_TFLA_S{i};
    StructArrayS0 = rac_TFLA_S0{i};
    for j = 1:length(StructArrayS)
        rac_TFLA_S_spectrum = StructArrayS(j).Data;
        rac_TFLA_S0_spectrum = StructArrayS0(j).Data;
        
        I_S = sum(rac_TFLA_S_spectrum(SB), 2); % integrals
        I_S = sum(I_S(:));

        I_S0 = sum(rac_TFLA_S0_spectrum(SB), 2); % integrals
        I_S0 = sum(I_S0(:));
        
        % Stack the spectra by adding an offset
        rac_TFLA_S_spectrum = rac_TFLA_S_spectrum + (i-1)*Offset;
        rac_TFLA_S0_spectrum = rac_TFLA_S0_spectrum + (i-1)*Offset;
        
        Integ_rac_S(end+1) = I_S; 
        Integ_rac_S0(end+1) = I_S0;

        % plot spectra
        figure(2)
        subplot(2,2,1)
        plot(rac_TFLA_S_spectrum(SB))
        title('rac-TFLA 19F-19F S CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e4 2.8e5])
        axis square
        
        hold on;

        subplot(2,2,2)
        plot(rac_TFLA_S0_spectrum(SB))
        title('rac-TFLA 19F-19F S0 CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e4 2.8e5])
        axis square

        hold on;
        
        % Plot noise
%         subplot(2,4,3)
%         plot(rac_TFLA_S_spectrum(NB))
%         title('rac-TFLA 19F-19F S NOISE')
%         xlabel('Freq / Hz')
%         ylabel('Relative intensity')
%         axis square
% 
%         hold on;
%         
%         subplot(2,4,4)
%         plot(rac_TFLA_S0_spectrum(NB))
%         title('rac-TFLA 19F-19F S0 NOISE')
%         xlabel('Freq / Hz')
%         ylabel('Relative intensity')
%         axis square
% 
%         hold on;
    end
end

for i = 1:length(S_TFLA_S)
    StructArrayS = S_TFLA_S{i};
    StructArrayS0 = S_TFLA_S0{i};
    for j = 1:length(StructArrayS)
        S_TFLA_S_spectrum = StructArrayS(j).Data;
        S_TFLA_S0_spectrum = StructArrayS0(j).Data;
        
        I_S = sum(S_TFLA_S_spectrum(SB), 2); % integrals
        I_S = sum(I_S(:));

        I_S0 = sum(S_TFLA_S0_spectrum(SB), 2); % integrals
        I_S0 = sum(I_S0(:));
        
        % Stack the spectra by adding an offset
        S_TFLA_S_spectrum = S_TFLA_S_spectrum + (i-1)*Offset;
        S_TFLA_S0_spectrum = S_TFLA_S0_spectrum + (i-1)*Offset;
        
        Integ_S_S(end+1) = I_S; 
        Integ_S_S0(end+1) = I_S0;

        % plot spectra
        figure(2)
        subplot(2,2,3)
        plot(S_TFLA_S_spectrum(SB))
        title('S-TFLA 19F-19F S CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e4 2.8e5])
        axis square
        
        hold on;

        subplot(2,2,4)
        plot(S_TFLA_S0_spectrum(SB))
        title('S-TFLA 19F-19F S0 CODEX')
        xlabel('Freq / Hz')
        ylabel('Relative intensity')
        xlim([0 2500])
        ylim([-1e4 2.8e5])
        axis square

        hold on;
        
%         % Plot noise
%         subplot(2,4,7)
%         plot(S_TFLA_S_spectrum(NB))
%         title('S-TFLA 19F-19F S NOISE')
%         xlabel('Freq / Hz')
%         ylabel('Relative intensity')
%         axis square
% 
%         hold on;
%         
%         subplot(2,4,8)
%         plot(S_TFLA_S0_spectrum(NB))
%         title('S-TFLA 19F-19F S0 NOISE')
%         xlabel('Freq / Hz')
%         ylabel('Relative intensity')
%         axis square
% 
%         hold on;
    end
end

set(gcf, 'PaperOrientation', 'landscape'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 29.7 21]); 
set(gcf, 'PaperSize', [29.7 21]); % A4 size in centimeters (landscape)
print(gcf, 'TFLA_CODEX_plots', '-dpdf', '-r300');

%% NORMALIZE S/S0 AND PLOT
n = length(Integ_rac_S);
if mod(n, 2) ~= 0
    Integ_rac_S = Integ_rac_S(1:end-1); % Remove the last element if odd
    Integ_rac_S0 = Integ_rac_S0(1:end-1); % Remove the last element if odd
end

q = length(Integ_S_S);
if mod(q, 2) ~= 0
    Integ_S_S = Integ_S_S(1:end-1); % Remove the last element if odd
    Integ_S_S0 = Integ_S_S0(1:end-1); % Remove the last element if odd
end

norm_CODEX_rac = Integ_rac_S ./ Integ_rac_S0;
norm_CODEX_S = Integ_S_S ./ Integ_S_S0;

% Plot the results
figure(3);
subplot(2,1,1)
plot(norm_CODEX_rac, '-o');
title('rac-TFLA CODEX Evolution')
xlabel('Mixing Time')
ylabel('S / S0')
grid on;
subplot(2,1,2)
plot(norm_CODEX_S, '-o');
title('S-TFLA CODEX Evolution')
xlabel('Mixing Time')
ylabel('S / S0')
grid on;

set(gcf, 'PaperOrientation', 'landscape'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 29.7 21]); 
set(gcf, 'PaperSize', [29.7 21]); % A4 size in centimeters (landscape)
print(gcf, 'TFLA_CODEX_decay_plots', '-dpdf', '-r300');

%% FIT MODELS
% see: https://pubs.acs.org/doi/full/10.1021/ja0603406



%% SAVEFIGS


%% FUNCTIONS
% see: https://pubs.acs.org/doi/full/10.1021/ja0603406
function RMSD = computeRMSD(I_sim, I_exp, N)
    RMSD = sqrt(sum((I_sim - I_exp)/N));
end

function k_ii = computeRateConstant(r_ij, theta_ij)

end

