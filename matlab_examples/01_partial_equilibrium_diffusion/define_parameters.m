function param = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
param.amin = -1;
param.amax = 20;
param.zmin = 0.8;
param.zmax = 1.2;

param.min = [param.amin, param.zmin];
param.max = [param.amax, param.zmax];


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

% Earnings parameters:
param.L = 1;
param.zmean = 1;
param.theta_z = 0.25;
param.sig_z = 0.01;


end