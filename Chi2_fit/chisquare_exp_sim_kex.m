% Chisquare fit to extract K_ex
% Author: Ettore Bartalucci, RWTH Aachen
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 04.04.24
% Project: CF3 self decoupling

% Load data from the .mat files
load('dataset1.mat');
load('dataset2.mat');

% Select one column of intensity2 randomly
random_index = randi(size(intensity2, 2));
selected_intensity2 = intensity2(:, random_index);

% Ensure that intensity1 and selected_intensity2 have the same size
min_length = min(length(intensity1), length(selected_intensity2));
intensity1 = intensity1(1:min_length);
selected_intensity2 = selected_intensity2(1:min_length);

% Define chi-square function
chi_square = @(data1, data2) sum((data1(:) - data2(:)).^2);

% Perform 2D chi-square minimization
chi_square_value = chi_square(intensity1, selected_intensity2);

disp(['Chi-square value between dataset 1 and selected dataset 2: ', num2str(chi_square_value)]);

% Plot 2D heat map of chi-square values
figure;
plot(omega(1:min_length), intensity1, 'b', omega(1:min_length), selected_intensity2, 'r');
xlabel('Frequency');
ylabel('Intensity');
title('Comparison of Intensity Profiles');
legend('Dataset 1', 'Randomly Selected Dataset 2');
