function [diff, G, ss] = stationary(x, G, p)

%% AGGREGATES

% Guess:
ss.r = x;

% Exogenous primitive:
ss.Z = p.Z;

% Macro block:
ss.w = (1-p.beta) * ss.Z.^(1/(1-p.beta)) .* ...
        (p.beta./(ss.r + p.delta)).^(p.beta/(1-p.beta));


%% VFI
for j = 1:p.Nz
for k = 1:p.Nr
    G.income_azj(:, j, k) = ss.r * G.a + ss.w(k) * p.zz(j);
    bound_azj(:, j, k) = p.u1(G.income_azj(:, j, k));
end
end
G.income = reshape(G.income_azj, [G.J, p.Nd]);
bound    = reshape(bound_azj,    [G.J, p.Nd]);

% State-constrained boundary conditions:
left_bound  = bound(G.grid(:, 1) == 0, :);
right_bound = bound(G.grid(:, 1) == 1, :);
for j = 1:p.Nd
    BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF';
    BC{1}.left.f  = @(points) left_bound(j)  * ones(size(points, 1), 1);
    BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
    G = gen_FD(G, BC, num2str(j));
end

% Initialize guess V0:
if ~isfield(G, 'V0'), G.V0 = p.u(G.income) / p.rho; end

% Solve VFI:
[ss.V, hjb] = HJB(G, p);
ss.c = hjb.c; ss.s = hjb.s; ss.u = hjb.u;
if(any(isnan(ss.V))), diff = NaN(1); return; end


%% OUTPUT VF AS NEXT GUESS
G.V0 = ss.V;


%% KOLMOGOROV FORWARD
ss.g = KF(hjb.Aa, hjb.Az, hjb.Api, G, p);


%% MARKET CLEARING

% Aggregation:
ss.l = sum(reshape(p.z.*ss.g, [p.Naz, p.Nr]) * G.da, 1);
ss.Q = sum(sum( G.a .* ss.g .* G.da));
ss.C = sum(sum( ss.c .* ss.g * G.da));
ss.S = sum(sum( ss.s .* ss.g .* G.da));

ss.mass = sum(reshape(ss.g, [p.Naz, p.Nr]), 1) * G.da;

% Remaining macro block:
ss.y = ss.l .* ss.Z.^(1/(1-p.beta)) .* (p.beta./(ss.r + p.delta)).^(p.beta/(1-p.beta));
ss.k = p.beta * ss.y ./ (ss.r + p.delta);

ss.K = sum(ss.k);
ss.I = p.delta*ss.K;

ss.Y = sum(ss.y);

% Some checks:
ss.w2 = (1-p.beta) .* ss.y ./ ss.l;

% Market clearing:
ss.excess_goods  = ss.Y - ss.C - ss.I;
ss.excess_wealth = ss.Q - ss.K;
ss.excess_saving = ss.S;

diff = ss.excess_wealth;

end




