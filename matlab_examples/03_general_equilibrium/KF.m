function g = KF(Aa, Az, G, p)

AT = (Aa + Az)';

% KF #1:
AT1 = AT;
b = zeros(p.discrete_types * G.J, 1);

i_fix = 1;
b(i_fix) = 0.1;
row = [zeros(1, i_fix-1), 1, zeros(1, p.discrete_types*G.J-i_fix)];
AT1(i_fix,:) = row;

gg = AT1 \ b;
g_sum = gg' * ones(p.discrete_types*G.J, 1) * G.da;
gg = gg ./ g_sum;

g1 = reshape(gg, [G.J, p.discrete_types]);


% KF #2:
g = zeros(G.J, p.discrete_types);
g(G.a == p.amin, :) = 1/p.discrete_types/G.da;

for n = 1:p.maxit_KF   
    B = 1/p.Delta_KF .* speye(p.discrete_types*G.J) - AT;
    b = g(:) / p.Delta_KF;

    g_new = B\b;

    diff = max(abs(g(:) - g_new));
    if diff < p.crit_KF, break; end
    g = reshape(g_new, [G.J, p.discrete_types]);
end
if n == p.maxit_KF, fprintf('KF did not converge. Remaining Gap: %.2d\n', diff); end

% Some tests:
mass = sum(g * G.da);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end
if max(max(abs(g1 - g))) > 1e-5, fprintf('Distributions g1 and g2 do not align!\n'); end


end