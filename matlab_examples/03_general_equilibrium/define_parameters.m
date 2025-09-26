function p = define_parameters(varargin)

%% GRID PARAMETERS

% Grid construction:
p.l = 8;
p.d = 1;

p.amin = 0;
p.amax = 100;

p.min = [p.amin];
p.max = [p.amax];


%% PDE TUNING PARAMETERS
p.Delta = 1000;
p.maxit = 100;
p.crit  = 1e-8;

p.Delta_KF = 1000;
p.maxit_KF = 100;
p.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS

% Household parameters:
p.rho = 0.02;
p.gamma = 2;
p.beta = 0.33;
p.delta = 0.025;

p.u     = @(x) x.^(1-p.gamma) / (1-p.gamma); 
p.u1    = @(x) x.^(-p.gamma);
p.u1inv = @(x) x.^(-1/p.gamma);

% Productivity: 
p.Z = 1;

% Earnings parameters:
p.z  = [0.8, 1.2];
p.la1 = 1/3;
p.la2 = 1/3;
p.L   = p.la2/(p.la1+p.la2) * p.z(1) + p.la1/(p.la1+p.la2) * p.z(2);
p.lambda = [-p.la1,  p.la1; p.la2, -p.la2]; % transition matrix for income
p.discrete_types = numel(p.z);

end