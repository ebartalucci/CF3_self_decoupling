% Fitting static T2 curves with integrals extracted manually from DMFIT

data = readtable('static_T2_dmfit_fits.xlsx');

delays = [1e-5,2.5e-5, 5e-5,7.5e-5, 1e-4, 1.5e-4, 2e-4, 2.5e-4, 3e-4,...
    3.5e-4, 4e-4, 4.5e-4, 5e-4, 6e-4,7.5e-4, 1e-3, 1.5e-3, 2e-3, 3e-3, 5e-3];

int_stfla1 = data.S_TFLA_p1;
int_stfla2 = data.S_TFLA_p2;
int_ractfla1 = data.rac_TFLA_p1;
int_ractfla2 = data.rac_TFLA_p2;

t = delays';
I = int_stfla2(:,1);

% Define the exponential decay function
expDecay = @(params, t) params(1) * exp(-t / params(2));

% Initial guesses for I0 and T2
initialGuess = [I(1), 0.001];  % [I0_guess, T2_guess]

% Fit the model using lsqcurvefit
options = optimset('Display', 'off');  % Suppress output display
params = lsqcurvefit(expDecay, initialGuess, t, I, [], [], options);

% Extract fitted parameters
I0 = params(1);
T2 = params(2);

% Display the results
disp(['I(0) = ', num2str(I0)]);
disp(['T2 = ', num2str(T2)]);

% Plot the data and the fitted curve
figure;
plot(t, I, 'bo', 'MarkerFaceColor', 'b');  % Original data points
hold on;
tFit = linspace(min(t), max(t), 100);  % Smooth curve for fit visualization
plot(tFit, expDecay(params, tFit), 'r-', 'LineWidth', 2);  % Fitted curve
xlabel('echo delay / s');
ylabel('intensity / a.u.');
title('T_2 Fit for broad S-TFLA');
legend('S-TFLA', 'Fit');