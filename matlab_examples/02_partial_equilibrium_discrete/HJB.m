function [V, c, s, u, A] = HJB(G, param)

V = G.V0;

% Generator for earnings process (2-state Markov chain):
Az = kron(param.lambda, speye(G.J));


for iter = 1:param.maxit
    
    % Solve for consumption policy function using upwind scheme:
    VaF = zeros(G.J, param.discrete_types);
    VaB = zeros(G.J, param.discrete_types);
    for j = 1:param.discrete_types
        VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
        VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
    end
    
    cF = param.u1inv(VaF);
    cB = param.u1inv(VaB);
    c0 = G.income;
    
    sF = G.income - cF;
    sB = G.income - cB;
    
    IF = (sF > 0);
    IB = (sB < 0) .* (IF == 0);
    I0 = (1-IF-IB);
    
    s = sF.*IF + sB.*IB;
    c = cF.*IF + cB.*IB + c0.*I0;
    u = param.u(c);

    % Compute generator for wealth evolution:
    Aa1 = FD_operator(G, s(:, 1), zeros(G.J, 1), 1, '1');
    Aa2 = FD_operator(G, s(:, 2), zeros(G.J, 1), 1, '2');
        
    % Solve linear system:
    A = blkdiag(Aa1, Aa2) + Az;
    B = (1/param.Delta + param.rho)*speye(G.J*param.discrete_types) - A;
    b = u(:) + V(:) / param.Delta;
    
    V_new = B\b;
    
    % Update value function:
    V_change = V_new - V(:);
    V = reshape(V_new, [G.J, param.discrete_types]);
    
    dist = max(max(abs(V_change)));
    if dist < param.crit; break; end    
    if mod(iter,1)==0, fprintf('VFI: %.i    Remaining Gap: %.2d\n',iter,dist); end
    
end

end
