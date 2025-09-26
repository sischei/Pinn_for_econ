function [diff, G, sim] = transition(PHI, G, shock, ss, p)

%% PREALLOCATE
sim.N = p.N; sim.t = p.t; sim.dt = p.dt;

[sim.V, sim.A, sim.Aa, sim.g, sim.c, sim.s] = deal(cell(p.N, 1));
[sim.Q, sim.C, sim.S, sim.r, sim.K, sim.dK, sim.K2] = deal(zeros(p.N, 1));


%% AGGREGATE TRANSITION PATH
X = basis_fun_irf([], reshape(PHI, [1, numel(PHI)]), p.H(1), p.H(2), ...
    p.bfun_type, sim.t, "get_function");

sim.r = X;

% Shock:
% if p.permanent_shock
%     sim.Z = p.ZT.*ones(p.N, 1);
% else
%     sim.Z = p.Z.*ones(p.N, 1) + [zeros(p.N, p.Nr-1), shock];
% end
sim.Z = p.Z.*ones(p.N, 1) + [zeros(p.N, p.Nr-1), shock]; % shock; %

% Macro block:
sim.w = (1-p.beta) .* sim.Z.^(1/(1-p.beta)) .* ...
        (p.beta./(sim.r + p.delta)).^(p.beta/(1-p.beta));


%% INITIALIZE FROM STEADY STATE
V = ss.V;
sim.g{1} = ss.g;


%% SOLVE HJB BACKWARDS
Az = kron(speye(p.Nr), kron(p.lambda, speye(G.J)));

for n = p.N:-1:1

    for j = 1:p.Nz
    for k = 1:p.Nr
        G.income_azj(:, j, k) = sim.r(n) * G.a + sim.w(n, k) * p.zz(j);
    end
    end
    G.income = reshape(G.income_azj, [p.Na, p.Nd]);
    
    % Consumption and savings:
    hjb = savings(V, G, p);

    % Migration:
    [pi, Api, V_gain] = migration(V, p);
    
    % Assemble FD operator matrix:
    Aa = [];
    for j = 1:p.Nd
        Aa = blkdiag(Aa, FD_operator(G, hjb.s(:, j), zeros(G.J, 1), 1, num2str(j)));
    end
    
    A = Aa + Az;
    
    B = (1/sim.dt(n) + p.rho)*speye(p.Ntot) - A - Api;
    b = hjb.u(:) + p.mu*V_gain + V(:) / sim.dt(n);
    
    % Solve linear system:
    V_new = B\b;
    V = reshape(V_new, [G.J, p.Nd]);
    
    if ~isreal(V), disp('Complex values detected!'); diff = NaN(1); return; end
    
    % Record data:
    sim.V{n} = V; sim.u{n} = hjb.u; sim.c{n} = hjb.c; sim.s{n} = hjb.s; sim.pi{n} = pi;
    sim.A{n} = A; sim.Aa{n} = Aa; sim.Api{n} = Api;
end


%% SOLVE KF FORWARDS
for n = 1:p.N
    
    AT = (sim.A{n} + sim.Api{n})';
    
    if p.implicit_g
        % Implicit:
        B = 1/sim.dt(n) * speye(p.Nd*G.J) - AT;
        b = sim.g{n}(:)/sim.dt(n);

        gg = B\b;
        sim.g{n+1} = reshape(gg, [G.J, p.Nd]);
    else
        % Explicit: (for explicit, N must be >6 times larger than T)
        gg = sim.g{n}(:) + sim.dt(n) * AT * sim.g{n}(:);
        sim.g{n+1} = reshape(gg, [G.J, p.Nd]);
    end
    sim.mass(n, :) = sum(sim.g{n} .* G.da);
    
    assert(all(all( sim.g{n+1} >= 0 )));   
    if abs( sum(sum(sim.g{n+1} .* G.da )) - 1) > 1e-8, error('Sim KF not mass-preserving.\n'); end 
    if abs( sum(sum(sim.g{n+1} .* G.da )) - sum(sum(sim.g{n} * G.da ))) > 1e-8
        fprintf('KF not preserving mass.\n'); end
end


%% AGGREGATION & MARKET CLEARING

% Aggregation:
for n = 1:p.N
    sim.l(n, :) = sum(reshape(p.z.*sim.g{n}, [p.Naz, p.Nr]) * G.da, 1);

    sim.Q(n) = sum(sum(G.a .* sim.g{n} .* G.da));
    sim.C(n) = sum(sum(sim.c{n} .* sim.g{n} .* G.da));
    sim.S(n) = sum(sum(sim.s{n} .* sim.g{n} .* G.da));
end

% Production:
sim.y = sim.l .* sim.Z.^(1/(1-p.beta)) .* (p.beta./(sim.r + p.delta)).^(p.beta/(1-p.beta));
sim.k = p.beta * sim.y ./ (sim.r + p.delta);

sim.K = sum(sim.k, 2);
sim.Y = sum(sim.y, 2);
for n = 1:p.N-1
    sim.dK(n) = (sim.K(n+1) - sim.K(n))/sim.dt(n);
end
sim.I = sim.dK + p.delta*sim.K;
sim.rk = sim.r + p.delta;

% Market clearing:
sim.excess_goods  = sim.Y - sim.C - sim.I;
sim.excess_wealth = sim.Q - sim.K;
sim.excess_saving = sim.S - sim.dK;


%% COLLOCATION POINTS
DIFF_Y = interp1(sim.t, sim.excess_goods,   p.nodes);
DIFF_K = interp1(sim.t, sim.excess_wealth,  p.nodes);
DIFF_S = interp1(sim.t, sim.excess_saving,  p.nodes);

diff = DIFF_Y';

end




