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
inv_mas_rates = 1./ mas_rates;

linfit_kex_no_apod_rac = polyval(polyfit(inv_mas_rates, kex_no_apod_tla_rac, 1), inv_mas_rates);

linfit_kex_no_apod_s = polyval(polyfit(inv_mas_rates, kex_no_apod_tla_s, 1), inv_mas_rates);

% Plot the values as inverted MAS
figure(1); hold on;

plot(inv_mas_rates, kex_no_apod_tla_rac, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
    '#A2142F', 'Color', '#A2142F')

plot(inv_mas_rates,kex_no_apod_tla_s, '-d', 'MarkerSize',4, 'LineWidth', 2, 'MarkerFaceColor',...
     '#EDB120', 'Color', '#EDB120')

plot(inv_mas_rates, linfit_kex_no_apod_rac, '--k', 'LineWidth', 1.2)

disp(['Equation for Rac is y = ' num2str(linfit_kex_no_apod_rac(1)) '*x + ' num2str(linfit_kex_no_apod_rac(2))])
SStot_rac = sum((kex_no_apod_tla_rac-mean(kex_no_apod_tla_rac)).^2);                    % Total Sum-Of-Squares
SSres_rac = sum((kex_no_apod_tla_rac-linfit_kex_no_apod_rac).^2);                       % Residual Sum-Of-Squares
Rsq_rac = 1-SSres_rac/SStot_rac;   

plot(inv_mas_rates, linfit_kex_no_apod_s, '--k', 'LineWidth', 1.2)

disp(['Equation for S is y = ' num2str(linfit_kex_no_apod_s(1)) '*x + ' num2str(linfit_kex_no_apod_s(2))])
SStot_s = sum((kex_no_apod_tla_s-mean(kex_no_apod_tla_s)).^2);                    % Total Sum-Of-Squares
SSres_s = sum((kex_no_apod_tla_s-linfit_kex_no_apod_s).^2);                       % Residual Sum-Of-Squares
Rsq_s = 1-SSres_s/SStot_s;                            % R^2

labels = string(mas_rates);
text(inv_mas_rates, kex_no_apod_tla_rac, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
axis('padded');

set(gca, 'YDir','reverse')
ylim([-10 800])

xlabel('1/\omega_{MAS} [s]')
ylabel('k_{ex} [s^{-1}]')

title('Variation of k_{ex} vs 1/\omega_{MAS}')

legend('rac-TFLA','({\itS})-TFLA', 'Location', 'southwest')

hold off;


