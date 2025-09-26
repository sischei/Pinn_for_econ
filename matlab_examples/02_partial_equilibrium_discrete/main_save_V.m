%------------------------------------------------------------------------%
% Model 02: solve discrete‐state HJB, dump plotting data, and plot lines  %
%------------------------------------------------------------------------%

clear
close all
clc

% ensure output folder exists
if ~exist('output','dir')
    mkdir('output')
end

diary ./output/output.log
diary on

addpath(genpath('../lib/'))
figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS
param = define_parameters();


%% INITIALIZE GRID
G = setup_grid(8, 0, param.min, param.max, ...
    'NamedDims', {1}, 'Names', {'a'});


%% SOLVE HJB

% Partial equilibrium prices
r = 0.01;
w = 1;

% Budget constraint: (J×2 matrix for two discrete states)
G.income = r * G.a + w .* param.zz;

% Boundary conditions for each discrete state
left_bound  = param.u1(G.income(1, :));
right_bound = param.u1(G.income(end, :));
for j = 1:param.discrete_types
    BC{1}.left.type  = 'VNB';
    BC{1}.right.type = 'VNF';
    BC{1}.left.f     = @(pts) left_bound(j)  * ones(size(pts,1),1);
    BC{1}.right.f    = @(pts) right_bound(j) * ones(size(pts,1),1);
    G = gen_FD(G, BC, num2str(j));
end

% Initial guess
G.V0 = param.u(G.income) / param.rho;

% Solve HJB
[V, c, s, u, A] = HJB(G, param);

fprintf('\nAlgorithm converged in %.2f seconds.\n', toc(run_time));


%% DUMP DATA FOR PLOTTING
% G.a is J×1, V,c,s are J×D where D = param.discrete_types (e.g. 2)
a_vec = G.a;
% Flatten each column
V1 = V(:,1);  V2 = V(:,2);
c1 = c(:,1);  c2 = c(:,2);
s1 = s(:,1);  s2 = s(:,2);

T = table(a_vec, V1, V2, c1, c2, s1, s2, ...
    'VariableNames', {'a', 'V1', 'V2', 'c1', 'c2', 's1', 's2'});
writetable(T, fullfile('output','plot_data_model02.csv'));
fprintf('Exported plotting data to output/plot_data_model02.csv\n');


%% PLOT FIGURES
figure; plot(G.a, V);   title('Value Function V(a)');      xlabel('a'); ylabel('V');   legend('z_1','z_2');
figure; plot(G.a, c);   title('Consumption c(a)');         xlabel('a'); ylabel('c');   legend('z_1','z_2');
figure; plot(G.a, s);   title('Savings s(a)');            xlabel('a'); ylabel('s');   legend('z_1','z_2');

diary off
