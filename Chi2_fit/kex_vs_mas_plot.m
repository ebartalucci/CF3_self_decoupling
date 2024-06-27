% Plot values of k_ex obtained from lstsq fit 
% vs MAS from the recorded spectrum 
% Author: Ettore Bartalucci, RWTH Aachen
% Project: CF3 self decoupling

%% Values
% K_ex values from the fits [s^-1]
kex_cos2_apod_tla_rac = [209 114.4 72.9 34.1 19.8 0 0];
kex_cos2_apod_tla_s = [716.5 534.5 354.2 190.2 85.6 12.9 4.4];

kex_no_apod_tla_rac = [188.6 102.5 72.2 48.8 34.5 0 0];
kex_no_apod_tla_s = [766.4 544.9 330.9 155.8 88.9 29.5 17.3];

mas_rates = [14 17.5 22 30 40 50 60];

%% Plots

% Plot the values as inverted MAS
figure(1); hold on;
inv_mas_rates = 1./ mas_rates;

plot(inv_mas_rates, kex_cos2_apod_tla_s, '-o', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#7E2F8E', 'Color', '#7E2F8E')
plot(inv_mas_rates, kex_cos2_apod_tla_rac, '-o', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
     '#0072BD', 'Color', '#0072BD')

plot(inv_mas_rates, kex_no_apod_tla_rac, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#A2142F', 'Color', '#A2142F')
plot(inv_mas_rates,kex_no_apod_tla_s, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
     '#EDB120', 'Color', '#EDB120')

set(gca, 'YDir','reverse')
ylim([-10 800])

xlabel('1/\omega_{MAS} [s]')
ylabel('k_{ex} [s^{-1}]')

title('SI figure k_{ex} vs 1/\omega_{MAS}')

legend('({\itS})-TFLA (cos^2 apodization)', 'rac-TFLA (cos^2 apodization)','({\itS})-TFLA', 'rac-TFLA', 'Location', 'southwest')

hold off;

% Plot the values as inverted MAS
figure(2); hold on;
inv_mas_rates = 1./ mas_rates;

plot(inv_mas_rates, kex_no_apod_tla_rac, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#A2142F', 'Color', '#A2142F')
plot(inv_mas_rates,kex_no_apod_tla_s, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
     '#EDB120', 'Color', '#EDB120')

labels = string(mas_rates);
dy = 0.1;
text(inv_mas_rates, kex_no_apod_tla_rac, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
axis('padded');

set(gca, 'YDir','reverse')
ylim([-10 800])

xlabel('1/\omega_{MAS} [s]')
ylabel('k_{ex} [s^{-1}]')

title('Variation of k_{ex} vs 1/\omega_{MAS}')

legend('({\itS})-TFLA', 'rac-TFLA', 'Location', 'southwest')

hold off;


