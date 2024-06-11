% Plot values of k_ex obtained from lstsq fit vs MAS from the recorded spectrum 
% Author: Ettore Bartalucci, RWTH Aachen
% Project: CF3 self decoupling

% K_ex values from the fits [s^-1]
kex_no_apod_tla_rac = [188.6 48.8 0];
kex_apod_tla_rac = [209 34.1 0];

kex_no_apod_tla_s = [766.4 155.8 17.3];
kex_apod_tla_s = [716.5 190.2 4.4];

mas_rates = [14 30 60];

% Plot the values
figure(1); clf; hold on;
plot(mas_rates, kex_apod_tla_s, '-o', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#7E2F8E', 'Color', '#7E2F8E')
plot(mas_rates, kex_no_apod_tla_s, '-o', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#0072BD', 'Color', '#0072BD')

plot(mas_rates, kex_apod_tla_rac, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#A2142F', 'Color', '#A2142F')
plot(mas_rates,kex_no_apod_tla_rac, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#EDB120', 'Color', '#EDB120')

set(gca, 'YDir','reverse')
ylim([-10 800])
xlim([10 62])

xlabel('MAS / kHz')
ylabel('k_{ex} / s^{-1}')

legend('({\itS})-TFLA (cos^2 apodization)', '({\itS})-TFLA', '({\itrac})-TFLA (cos^2 apodization)', '({\itrac})-TFLA', 'Location', 'southeast')

hold off;
