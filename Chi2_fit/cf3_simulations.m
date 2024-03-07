% Chi^2 Fit for the CF3 lineshape project
% Author: Ettore Bartalucci
% Aachen, 07.03.24
% Function for simulating CF3 lineshape
% Adapted from Mathematica script of Prof. Matthias Ernst

%% Initialize variables
k_ex = 1/1000;
T_2 = 3/100;
J_cf = 280; % Hz
I = 1/2;


%% Get FID

L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex; k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit

L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
      3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
      0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
      0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];

% Compute and simplify the matrix exponential
syms t;
U1 = simplify(expm(L1 * abs(t))); 
U2 = simplify(expm(L2 * abs(t)));

% Compute signal1
signal1 = sum(U1 * [2; 2]);

% Plotting
t_values = linspace(0, 0.1, 1000); % Generate 1000 time values from 0 to 0.1
signal_values = double(subs(signal1, [t, T2, k0], [t_values.', T2, k0])); % Substitute t, T2, and k0
plot(t_values, signal_values);
xlabel('t');
ylabel('Signal1');
title('Plot of Signal1');



