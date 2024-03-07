% Chi^2 Fit for the CF3 lineshape project
% Author: Ettore Bartalucci
% Aachen, 07.03.24





% Generate example data
xData = rand(100,1);
yData = rand(100,1);
zData = 2 + 3*xData + 4*yData + randn(100,1);

% Initial guess for parameters
initialGuess = [1, 1, 1];

% Perform the fit
[fitParams, chiSquare] = chiSquareFit2D(xData, yData, zData, initialGuess);

% Display results
disp('Fit parameters:');
disp(fitParams);
disp(['Chi-square value: ', num2str(chiSquare)]);