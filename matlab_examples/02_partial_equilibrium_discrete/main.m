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
param = define_parameters();


%% INITIALIZE GRIDS
G = setup_grid(8, 0, param.min, param.max, 'NamedDims', {1}, 'Names', {'a'});


%% SOLVE HJB

% Partial equilibrium:
r = 0.01;
w = 1;

% Budget constraint:
G.income = r * G.a + w .* param.zz;

% Boundary conditions:
left_bound  = param.u1(G.income(G.grid(:, 1) == 0, :));
right_bound = param.u1(G.income(G.grid(:, 1) == 1, :));
for j = 1:param.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

% Initialize guess for value function:
G.V0 = param.u(G.income) / param.rho;

% Solve HJB:
[V, c, s, u, A] = HJB(G, param);


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

figure; plot(G.a, V);
figure; plot(G.a, c);
figure; plot(G.a, s);

diary off




