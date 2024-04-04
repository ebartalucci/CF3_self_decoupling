% Test work with NMR data in matlab
% Author: Ettore Bartalucci, RWTH Aachen
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 04.04.24
% Project: CF3 self decoupling

%% Open spectrum using RBNMR function and plot to see if everything alright
NMR_data = rbnmr('D:\PhD\Publications\CF3_MAS_SD_Jcoupling\codes\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');

%% Select region in the CF3 region of the spectrum (115ppm-135ppm)
w = NMR_data.XAxis;
spectrum = NMR_data.Data;
w_cf3 = w(3400:3700);
spectrum_cf3 = spectrum(3400:3700);

figure(1); clf;
subplot(2,1,1)
plot(spectrum)
title('Spectrum in number of points')
xlabel('N of points')

subplot(2,1,2)
plot(w, spectrum)
set(gca, 'XDir', 'reverse');
title('Spectrum in ppm')
xlabel('13C / ppm')

figure(2); clf;
subplot(2,1,1)
plot(spectrum_cf3)
title('CF3 region in number of points')
xlabel('N of points')

subplot(2,1,2)
plot(w_cf3, spectrum_cf3)
set(gca, 'XDir', 'reverse');
title('CF3 region in ppm')
xlabel('13C / ppm')


save('exp_cf3_data_TLA_S_14khz_exp_200_161023.xlsx', 'w_cf3', "spectrum_cf3")


