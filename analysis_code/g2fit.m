%%
% Function for fitting g(2)(t) according to a three-level system
%%

function [fitresult, gof] = g2fit(tau, g2y, param_init)
if nargin < 3
    param_init = [0.8 1 15e-9 250e-9 0];
end

[xData, yData] = prepareCurveData( tau, g2y );

% Set up fittype and options.
ft = fittype( '(1-a*(b*exp(-abs(tau-e)/c)+(1-b)*exp(-abs(tau-e)/d)))', 'independent', 'tau', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'notify';
opts.Lower = [0 0 0 0 -Inf];
opts.StartPoint = param_init;
opts.Upper = [1 Inf Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end