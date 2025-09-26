function g = KF(Aa, Az, Api, G, p)

AT = (Aa + Az + Api)';

% KF #1:
% AT1 = AT;
% b = zeros(p.Ntot, 1);
% 
% i_fix = 1;
% b(i_fix) = 0.1;
% row = [zeros(1, i_fix-1), 1, zeros(1, p.Ntot-i_fix)];
% AT1(i_fix, :) = row;
% 
% gg = AT1 \ b;
% g_sum = gg' * ones(p.Ntot, 1) * G.dx;
% gg = gg ./ g_sum;
% 
% g1 = [gg(1:G.J), gg(G.J+1:2*G.J)];

% KF #2:
g = zeros(G.J, p.Nd);
g(G.a == p.amin, :) = 1/p.Nd/G.da;

for n = 1:p.maxit_KF   
    B = 1/p.Delta_KF .* speye(p.Ntot) - AT;
    b = g(:)/p.Delta_KF;
    
    g_new = B\b;
    
    diff = max(abs( g(:) - g_new ));
    if diff < p.crit_KF, break; end
    g = reshape(g_new, [G.J, p.Nd]);
end
if n == p.maxit_KF, fprintf('KF did not converge. Remaining Gap: %.2d\n', diff); end

% Some tests:
mass = sum(g * G.da);
if abs(sum(mass)-1) > 1e-5, fprintf('Distribution not normalized!\n'); end
% if max(max(abs(g1 - g))) > 1e-5, fprintf('Distributions g1 and g2 do not align!\n'); end

end