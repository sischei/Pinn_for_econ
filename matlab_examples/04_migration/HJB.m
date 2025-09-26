function [V, hjb] = HJB(G, p)

V = G.V0;

% Exogenous operators:
Az = kron(speye(p.Nr), kron(p.lambda, speye(G.J)));

for iter = 1:p.maxit
    
    % Consumption and savings:
    hjb = savings(V, G, p);
    if any(any(isnan(hjb.c))), V = NaN(1); return; end
    
    % Migration:
    [pi, Api, V_gain] = migration(V, p);

    % Assemble FD operator matrix:
    Aa = [];
    for j = 1:p.Nd
        Aa = blkdiag(Aa, FD_operator(G, hjb.s(:, j), zeros(G.J, 1), 1, num2str(j)));
    end
    
    A = Aa + Az;
    
    B = (1/p.Delta + p.rho)*speye(p.Ntot) - A - Api;
    b = hjb.u(:) + V(:) / p.Delta + p.mu*V_gain;
    
    % Solve linear system:
    V_new = B\b;
    
    V_change = V_new - V(:);
    V = reshape(V_new, [G.J, p.Nd]);
    
    dist = max(max(abs(V_change)));
    if dist < p.crit; break; end
    
    % if mod(iter, 1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end
    if ~isreal(V), fprintf('Complex values in VFI: terminating process.\n'); V = NaN(1); return; end

end

hjb.A = A; hjb.Aa = Aa; hjb.Az = Az; hjb.Api = Api; hjb.pi = pi;

% Check for convergence:
if iter == p.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
