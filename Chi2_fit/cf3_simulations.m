% Chi^2 Fit for the CF3 lineshape project
% Author: Ettore Bartalucci
% Aachen, 07.03.24
% Function for simulating CF3 lineshape
% Adapted from Mathematica script of Prof. Matthias Ernst

tic;

%% Initialize variables
k_ex = 1/1000;
T_2 = 3/100;
J_cf = 280; 

%% Simulate FID
% Liouvillian superoperators. The matrices contain the dynamics of the spin system
L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf), k_ex;
     k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf)]; %1i is the imaginary unit

L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf), 3*k_ex, 0, 0;
      3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf), 4*k_ex, 0;
      0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf), 3*k_ex;
      0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf)];

% Compute and simplify the matrix exponential to get the time evolution
sym t
U1 = expm(L1 * abs(t));
U1 = simplify(U1);

% U2 = expm(L2 * abs(t));
% U2 = simplify(U2);

% Get the signals
signal1 = sum(U1 * [2; 2]);
% signal2 = sum(U2 * [1; 1; 1; 1]);

% Plotting
t_vals = linspace(0, 0.1, 50); % Generate 50 time values from 0 to 0.1
signal1_values = double(subs(signal1, [t, T_2, k_ex], [t_value, T_2, k_ex])); % Substitute t, T2, and k0
% signal2_values = double(subs(signal2, [t, T_2, k_ex], [t_value, T_2, k_ex]));

figure(1);clf;
plot(t_value, signal1_values);
xlabel('t');
ylabel('Signal1');
title('Plot of Signal1');

% figure(2);clf;
% plot(t_value, signal2_values);
% xlabel('t');
% ylabel('Signal2');
% title('Plot of Signal2');
% 





runtime = toc;
disp(['Elapsed time: ', num2str(runtime), ' seconds']);





