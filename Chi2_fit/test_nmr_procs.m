%% Experimental section
% Load the experimental spectra
% Optical pure @ 14kHz, 30kHz, and 60kHz
tla_s_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_14khz_exp_200_161023\pdata\1');
tla_s_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_30khz_exp_14_121023\pdata\1');
tla_s_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_S_60khz_exp_15_121023\pdata\1');
% Racemic @ 14kHz, 30kHz, and 60kHz
tla_rac_14khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_14khz_exp_10_221123\pdata\1');
tla_rac_30khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_30khz_exp_10_211123\pdata\1');
tla_rac_60khz = rbnmr('C:\Users\ettor\Desktop\Programming_for_NMR_spectroscopists\CF3_self_decoupling\Chi2_fit\spectra\TLA_rac_60khz_exp_13_211123\pdata\1');

plot(tla_rac_14khz.Data)


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


