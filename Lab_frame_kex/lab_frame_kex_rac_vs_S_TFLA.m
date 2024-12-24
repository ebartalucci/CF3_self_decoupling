clf;
%% AVOID RELOADING
SkipVars = {'FOLDER'};

for i = 1:length(SkipVars)
    varName = SkipVars{i};
    if evalin('base', ['exist(''', varName, ''', ''var'')'])
        disp([varName ' already exists in the workspace, skipping...']);
    else
        disp([varName ' is not in the workspace, loading...']);
        FOLDER = rbnmr('D:\PhD\Data\NMR\ETBA\Bolm\CF3_coupling\061224_3p2mm_19F_SD_curves_TFLA_rac_S_new_TWC_EB');
    end
end

Spectra_S = FOLDER(56:71);

Spectra_rac = FOLDER(75:90);

SB_CF3 = 57000:65000;
NB = 90000:100000;

vdlist = [1e-6 2e-6 5e-6 10e-6 25e-6 50e-6 75e-6 100e-6 200e-6 350e-6 500e-6...
    700e-6 850e-6 1000e-6 1500e-6 2000e-6];

Offset = 100000; % for stacked plots

Integ_CF3_S = [18411732 21765953 21078448 18797899 19090434 16996804 ...
    18613949 17114785 16427294 14655169 13385137 7538818 5899239 ...
    3305658 551432 -1372532]; 

Integ_CF3_S = Integ_CF3_S./Integ_CF3_S(1);


Integ_CF3_rac = [21529573 21877111 18071405 17985929 16616788 17728855 ...
    20486527 21516983 20036695 12688804 15869424 8532690 5630354 4498845 ...
    -2098659 1629241]; 

Integ_CF3_rac = Integ_CF3_rac./Integ_CF3_rac(1);

%% GET SPECTRA
for i = 1:length(Spectra_S)
    StructArrayS = Spectra_S{i};
    StructArrayrac = Spectra_rac{i};

    S_spectrum = StructArrayS.Data;
    rac_spectrum = StructArrayrac.Data;

    % Stack the spectra by adding an offset
%     I_S_CF3 = sum(S_spectrum, 2); % integrals
%     I_S_CF3 = sum(I_S_CF3(:));
%     Integ_CF3_S(end+1) = I_S_CF3;

%     Integ_CF3_S = Integ_CF3_S./Integ_CF3_S(1);
% 
%     I_rac_CF3 = sum(rac_spectrum, 2); % integrals
%     I_rac_CF3 = sum(I_rac_CF3(:));
%     Integ_CF3_rac(end+1) = I_rac_CF3;
% 
%     Integ_CF3_rac = Integ_CF3_rac./Integ_CF3_rac(1);

    S_spectrum = S_spectrum(SB_CF3) + (i-1)*Offset;
    rac_spectrum = rac_spectrum(SB_CF3) + (i-1)*Offset;

    % plot spectra
    figure(1);
    subplot(1,3,1)
    plot(S_spectrum)
    title('(S)-TFLA Lab frame 19F SD')
    xlabel('Freq / Hz')
    ylabel('Relative intensity')
    xlim([2000 7000])
    
    hold on;

    subplot(1,3,2)
    plot(rac_spectrum)
    title('rac-TFLA Lab frame 19F SD')
    xlabel('Freq / Hz')
    ylabel('Relative intensity')
    xlim([2000 7000])
    
    hold on;

end

subplot(1,3,3)
plot(vdlist, Integ_CF3_S, '-o', Color='k')  
hold on;
plot(vdlist, Integ_CF3_rac, '-d', Color='r')  
legend('(S)-TFLA', 'rac-TFLA')
xlabel('\tau_{mix} / s')
ylabel('Normalized Intensity / a.u.')
axis square


set(gcf, 'PaperOrientation', 'landscape'); 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 29.7 21]); 
set(gcf, 'PaperSize', [29.7 21]); % A4 size in centimeters (landscape)
print(gcf, 'plots', '-dpdf', '-r300');