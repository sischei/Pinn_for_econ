function [V, policies] = HJB(G, p)

V = G.V0;

% Exogenous Operators:
Az = kron(p.lambda, speye(G.J));

for iter = 1:p.maxit

    % Consumption and savings:
    policies = savings(V, G, p);
    if any(any(isnan(policies.c))), V = NaN(1); return; end

    % Assemble FD operator matrix:
    Aa = [];
    for j = 1:p.discrete_types
        Aa = blkdiag(Aa, FD_operator(G, policies.s(:, j), zeros(G.J, 1), 1, num2str(j)));
    end

    A = Aa + Az;

    B = (1/p.Delta + p.rho)*speye(p.discrete_types*G.J) - A;
    b = policies.u(:) + V(:) / p.Delta;

    % Solve linear system:
    V_new = B\b;

    V_change = V_new - V(:);
    V = reshape(V_new, [G.J, p.discrete_types]);
    
    dist = max(max(abs(V_change)));
    if dist < p.crit; break; end
    
    % if mod(iter, 1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end
    if ~isreal(V), fprintf('Complex values in VFI: terminating process.\n'); V = NaN(1); return; end

end

policies.A = A; policies.Aa = Aa; policies.Az = Az;

% Check for convergence:
if iter == p.maxit, fprintf('VFI did not converge. Remaining Gap: %.2d\n', iter, dist); V = NaN(1); return; end

end
