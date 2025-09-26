function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.amin = -1;
param.amax = 20;

param.min = [param.amin];
param.max = [param.amax];


%% PDE TUNING PARAMETERS
param.Delta = 1000;
param.maxit = 100;
param.crit  = 1e-8;


%% ECONOMIC PARAMETERS

% Household parameters:
param.rho = 0.02;
param.gamma = 2;

param.u     = @(x) x.^(1-param.gamma) / (1-param.gamma); 
param.u1    = @(x) x.^(-param.gamma);
param.u1inv = @(x) x.^(-1/param.gamma);

% Earnings process:
param.zz = [0.8, 1.2];
param.lambda1 = 1/3;
param.lambda2 = 1/3;
param.lambda = [-param.lambda1,  param.lambda1; param.lambda2, -param.lambda2];
param.discrete_types = numel(param.zz);

end