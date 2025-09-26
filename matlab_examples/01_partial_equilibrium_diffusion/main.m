%------------------------------------------------------------------------%
% 
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
G = setup_grid(0, [8, 4], param.min, param.max, 'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);


%% COMPUTE VALUE FUNCTION

% Partial equilibrium:
r = 0.01; 
w = 1;

% Household budget constraint:
G.income = r * G.a + w .* G.z;

% Boundary conditions:
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);

BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
BC{1}.left.f  = @(points) sparse_project(left_bound,  points, G);
BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
BC{2}.left.type = '0'; BC{2}.right.type = '0';
G = gen_FD(G, BC);

% Initialize guess V0:
G.V0 = param.u(G.income) / param.rho;

% Solve VFI:
[V, c, s, u, A] = HJB(G, param);


%% OUTPUT
run_time = toc(run_time); fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

aa = unique(G.a); zz = unique(G.z); Ja = numel(aa); Jz = numel(zz);
figure; mesh(aa, zz, reshape(V, [Ja, Jz])');
figure; mesh(aa, zz, reshape(c, [Ja, Jz])');
figure; mesh(aa, zz, reshape(s, [Ja, Jz])');

diary off




