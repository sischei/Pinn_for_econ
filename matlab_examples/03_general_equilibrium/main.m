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


%% SOLVE STATIONARY EQUILIBRIUM
K0 = 30; X0 = K0;

% Get better guess for value function:
[diff0, G, ~] = stationary(X0, G, p);

% Solve for steady state prices:
options = optimset('Display', 'iter', 'UseParallel', false, 'TolX', 1e-14);
X = fsolve(@(x) stationary(x, G, p), X0, options);

% Solve with correct prices:
[~, G, ss] = stationary(X, G, p);

fprintf('Stationary Equilibrium: (r = %.4f, K = %.2f),  markets(S = %.2d,  Y-C-I = %.2d, Kgap = %.2d) \n\n', ...
    ss.r, ss.K, ss.S, ss.excess_supply, ss.excess_capital);



%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');    
figure('visible', 'on'); hold on;
l1=scatter(G.a, ss.V(:, 1)); 
l2=scatter(G.a, ss.V(:, 2)); 
hold off; xlabel('Capital');
legend([l1,l2], {'$V^U(k)$', '$V^E(k)$'}, 'Interpreter', 'Latex', 'box', 'off', 'Location', 'SouthEast');

diary off


