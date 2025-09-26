%------------------------------------------------------------------------%
%                                                                        %
%                      MAIN SCRIPT - HJB SOLVER                           %
%                                                                        %
%------------------------------------------------------------------------%

clear
close all
clc

% --- Setup ---
diary ./output/output.log
diary on

addpath(genpath('../lib/'))
% figure_format; % Optional custom formatting

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS
% NOTE: Assumes a function 'define_parameters' exists and returns a struct.
param = define_parameters();


%% INITIALIZE GRIDS
% NOTE: Assumes a function 'setup_grid' exists.
G = setup_grid(0, [8, 4], param.min, param.max, 'NamedDims', {1, 2}, 'Names', {'a', 'z'}, 'DxxDims', 2);


%% COMPUTE VALUE FUNCTION

% Partial equilibrium prices
r = 0.01;
w = 1;

% Household budget constraint
G.income = r * G.a + w .* G.z;

% Boundary conditions
% NOTE: Assumes param.u1 is defined in 'define_parameters'.
left_bound  = param.u1(G.income);
right_bound = param.u1(G.income);

% NOTE: Assumes functions 'sparse_project' and 'gen_FD' exist.
BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
BC{1}.left.f  = @(points) sparse_project(left_bound,  points, G);
BC{1}.right.f = @(points) sparse_project(right_bound, points, G);
BC{2}.left.type = '0'; BC{2}.right.type = '0';
G = gen_FD(G, BC);

% Initialize guess V0
% NOTE: Assumes param.u is defined in 'define_parameters'.
G.V0 = param.u(G.income) / param.rho;

% Solve VFI using a Hamilton-Jacobi-Bellman solver
% NOTE: Assumes a function 'HJB' exists.
[V, c, s, u, A] = HJB(G, param);


%% --- OUTPUT AND PLOTTING ---
run_time = toc(run_time);
fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

fprintf('\nPlotting Figures...\n');

% Get unique grid points for plotting axes
aa = unique(G.a);
zz = unique(G.z);
Ja = numel(aa);
Jz = numel(zz);

% Generate original Matlab figures
figure; mesh(aa, zz, reshape(V, [Ja, Jz])'); title('Value Function (V)');
figure; mesh(aa, zz, reshape(c, [Ja, Jz])'); title('Consumption Policy (c)');
figure; mesh(aa, zz, reshape(s, [Ja, Jz])'); title('Savings Policy (s)');


%------------------------------------------------------------------------%
%                  NEW: EXPORT ALL DATA TO ONE CSV FILE                  %
%------------------------------------------------------------------------%

fprintf('\nSaving all policy function data to CSV file...\n');

% Create a meshgrid of the 'a' and 'z' coordinates
[A_mesh, Z_mesh] = meshgrid(aa, zz);

% Reshape ALL policy vectors into matrices that match the grid dimensions.
% The transpose (') is crucial to align with how mesh() interprets data.
V_matrix = reshape(V, [Ja, Jz])';
C_matrix = reshape(c, [Ja, Jz])';
S_matrix = reshape(s, [Ja, Jz])';

% Flatten the matrices into single columns for the table
a_col = A_mesh(:);
z_col = Z_mesh(:);
v_col = V_matrix(:);
c_col = C_matrix(:);
s_col = S_matrix(:);

% Create a table with columns for all exported variables
T = table(a_col, z_col, v_col, c_col, s_col, 'VariableNames', {'a', 'z', 'V', 'c', 's'});

% Write the table to a single CSV file
writetable(T, 'HJB_data.csv');

fprintf('Successfully saved data to HJB_data.csv\n');

diary off
