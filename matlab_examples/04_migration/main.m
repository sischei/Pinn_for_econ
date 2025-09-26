%------------------------------------------------------------------------%
% 
%------------------------------------------------------------------------%

clear
close all
clc

diary ./output/output.log
diary on

addpath(genpath('../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS

p = define_parameters();


%% INITIALIZE GRIDS
G = setup_grid(p.l, 0, p.min, p.max, 'NamedDims', {1}, 'Names', {'a'});

% Grid points:
p.Na   = G.J;                   % assets
p.Nz   = numel(p.zz);           % income
p.Nd   = p.Nz * p.Nr;           % total discrete types
p.Naz  = G.J * p.Nz;            % grid points \ regions
p.Ntot = p.Na * p.Nz * p.Nr;    % all: assets and income


%% COMPUTE STATIONARY EQUILIBRIUM
r0 = 0.02; X0 = r0;

% Get better guess for value function:
[~, G, ~] = stationary(X0, G, p);

% Solve for steady state prices:
options = optimset('Display', 'iter', 'UseParallel', false, 'TolX', 1e-12);
X = fsolve(@(x) stationary(x, G, p), X0, options);

% Solve with correct prices:
[~, G, ss] = stationary(X, G, p);

fprintf('Stationary Equilibrium: (r=%.4f  Y=%.4f)  markets(S=%.2d  Y-C-I=%.2d  A-K=%.2d) \n\n', ...
    ss.r, ss.Y, ss.S, max(abs(ss.excess_goods)), ss.excess_wealth);
fprintf('Hash #1: %.12f\n', ss.r);
fprintf('Hash #2: %.12f\n', ss.K);



%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

figure; plot(G.a, ss.g);

diary off
