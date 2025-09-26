function p = define_parameters(varargin)

%% GRID PARAMETERS
p.l = 7;
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
p.maxit_KF = 300;
p.crit_KF  = 1e-8;


%% ECONOMIC PARAMETERS

% Migration:
p.Nr    = 5;   % number of regions
p.mu    = 1; % migration rate
p.theta = 0.5; % migration elasticity

temp     = 1:p.Nr; 
p.kappa  = -(temp - temp').^2; % migration costs
p.ekappa = exp(p.kappa);
p.etk    = exp(p.theta*p.kappa);

% Productivity shock:
p.Z = linspace(0.8, 1, p.Nr); % linspace(1, 1, p.Nr); % location productivities

% Preferences:
p.rho = 0.02;
p.gamma = 2;

p.u     = @(x) x.^(1-p.gamma) / (1-p.gamma); 
p.u1    = @(x) x.^(-p.gamma);
p.u1inv = @(x) x.^(-1/p.gamma);

% Production:
p.beta = 0.33;

% Earnings parameters:
p.zz  = [0.8, 1.2];
p.la1 = 1/3;
p.la2 = 1/3;
p.L   = p.la2/(p.la1+p.la2) * p.zz(1) + p.la1/(p.la1+p.la2) * p.zz(2);

p.z = repmat(p.zz, [1, p.Nr]);

p.lambda = [-p.la1,  p.la1; p.la2, -p.la2]; % transition matrix for income

p.Nz = numel(p.zz); % number of income states

% Capital construction:
% p.delta = 0.025 * ones(1, p.Nr);
p.delta = 0.025;

% Number of grid points:
p.Nd = p.Nr * p.Nz ; % discrete types: regions and income


%% TRANSITION DYNAMICS PARAMETERS
p.time_grid_adjustment = 1;
p.T = 100; 
p.N = 120;
p.t = linspace(0, p.T, p.N)';
%{
    Here is what works: T=100, N=200, "nodal".
    For "cheb", here is what works: T = 100; NN = 100 / 200; H >= 25; it seems "cheb" does not like sparser time grids
    In some cases I need to tune H a little for "cheb" (what usually works: 20 - 25)
%}

p.bfun_type = "cheb"; 
p.cheb_H = 25;

if p.time_grid_adjustment == 1
    if p.N / p.T >= 2
        adjustment = @(x) x; 
    elseif p.N / p.T >= 1
        adjustment = @(x) (exp(x/p.T)-1) * p.T / (exp(1)-1);
    elseif p.N / p.T > 0.8
        adjustment = @(x) x.^2 / p.T^1;
    else 
        adjustment = @(x) x.^3 / p.T^2;
    end
    p.t = adjustment(p.t);
end
p.dt = diff(p.t); p.dt(p.N) = p.dt(p.N-1);

p.H(1) = p.N; if p.bfun_type == "cheb", p.H(1) = p.cheb_H; end
p.H(2) = 1; % # of time series to guess 

if p.N <= 6*p.T, p.implicit_g = 1; else p.implicit_g = 0; end


%% VARIABLE INPUTS

% Parse inputs:
q = inputParser;
q.CaseSensitive = true;
for f = fieldnames(p)'
    q.addParameter(f{:}, p.(f{:}));
end
parse(q, varargin{:});
p = q.Results;

% Simulation:
p.reso_sim_KF = 7 - floor(p.N / p.T);
p.t = linspace(0, p.T, p.N)';

if p.time_grid_adjustment == 1
    if p.N / p.T >= 2
        adjustment = @(x) x; 
    elseif p.N / p.T >= 1
        adjustment = @(x) (exp(x/p.T)-1) * p.T / (exp(1)-1);
    elseif p.N / p.T > 0.8
        adjustment = @(x) x.^2 / p.T^1;
    else 
        adjustment = @(x) x.^3 / p.T^2;
    end
    p.t = adjustment(p.t);
end
p.dt = diff(p.t); p.dt(p.N) = p.dt(p.N-1);
if p.bfun_type == "cheb"
    p.H(1) = p.cheb_H; 
elseif p.bfun_type == "nodal"
    p.H(1) = p.N; 
end

if p.N < 6*p.T, p.implicit_g = 1; else p.implicit_g = 0; end

end