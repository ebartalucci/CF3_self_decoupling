% Chisquare fit function
% Author: Ettore Bartalucci, RWTH Aachen
% Support and debug with Chatgpt
% First draft: Aachen, 07.03.24
% Last update: Aachen, 04.04.24
% Project: CF3 self decoupling

function [fitParams, chiSquare] = chiSquareFit2D(xData, yData, zData, initialGuess)
    % xData, yData: arrays containing the x and y coordinates of the data points
    % zData: array containing the measured z values corresponding to the data points
    % initialGuess: initial guess for the parameters of the fit
    
    % Define the model function
    model = @(params, x, y) params(1) + params(2)*x + params(3)*y;
    
    % Define the chi-square function
    chiSquareFunc = @(params) sum((model(params, xData, yData) - zData).^2);
    
    % Use fminsearch to minimize the chi-square function
    fitParams = fminsearch(chiSquareFunc, initialGuess);
    
    % Calculate the chi-square value
    chiSquare = chiSquareFunc(fitParams);
    
    % Generate a grid of x and y values for plotting
    xValues = linspace(min(xData), max(xData), 100);
    yValues = linspace(min(yData), max(yData), 100);
    [X, Y] = meshgrid(xValues, yValues);
    
    % Calculate the fitted z values for the grid
    Z_fit = model(fitParams, X, Y);
    
    % Plot the data points
    figure(100);clf;
    scatter(xData, yData, 30, zData, 'filled');
    hold on;
    
    % Plot the contour plot of the fitted surface
    contour(X, Y, Z_fit, 2000, 'LineWidth', 2);
    colorbar;
    
    % Set labels and title
    xlabel('x');
    ylabel('y');
    title('2D Chi-Square Fit');
    
    % Display results
    disp('Fit parameters:');
    disp(fitParams);
    disp(['Chi-square value: ', num2str(chiSquare)]);
    
    hold off;
end
