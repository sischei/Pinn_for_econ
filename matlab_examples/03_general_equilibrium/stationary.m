function [diff, G, ss] = stationary(x, G, p)

%% AGGREGATES
K = x;

L  = p.L;
Y  = p.Z * K^p.beta * L^(1-p.beta);
rk = p.Z * p.beta * Y / K;
w  = p.Z * (1-p.beta) * Y / L;
r  = rk - p.delta;


%% VFI
G.r = r;
G.w = w;
G.income = r * G.a + w * p.z;

% State-constrained boundary conditions:
left_bound  = p.u1(G.income(G.grid(:, 1) == 0, :));
right_bound = p.u1(G.income(G.grid(:, 1) == 1, :));

for j = 1:p.discrete_types
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j)  * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G,'V0'), G.V0 = p.u(G.income) / p.rho; end

% Solve HJB:
[V, policies] = HJB(G, p);


%% OUTPUT VF AS NEXT GUESS
G.V0 = V;


%% KOLMOGOROV FORWARD
g = KF(policies.Aa, policies.Az, G, p);


%% MARKET CLEARING
KH = sum(sum(G.a .* g .* G.da));
LH = sum(sum(p.z .* g .* G.da));
C  = sum(sum(policies.c .* g .* G.da));
S  = sum(sum(policies.s .* g .* G.da));

excess_supply  = Y - C - p.delta*KH;
excess_capital = K - KH;

diff = excess_capital;

ss.V = V; ss.g = g; ss.c = policies.c; ss.s = policies.s;
ss.K = K; ss.L = L; ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w; 
ss.excess_supply = excess_supply; ss.excess_capital = excess_capital;

end




