function [V, c, s, u, A] = HJB(G, param)

V = G.V0;

% Generator for earnings process (OU diffusion):
Az = FD_operator(G, param.theta_z * (param.zmean - G.z), param.sig_z*ones(G.J,1), 2);

for iter = 1:param.maxit

    % Consumption policy function via upwind scheme:
    VaF = deriv_sparse(G, V, 1, 'D1F');
    VaB = deriv_sparse(G, V, 1, 'D1B');
    
    cF = param.u1inv(VaF);
    cB = param.u1inv(VaB);
    c0 = G.income;
    
    sF = G.income - cF;
    sB = G.income - cB;
    
    IF = (sF > 0) .* (G.grid(:, 1) < 1);  
    IB = (sB < 0) .* (IF == 0) .* (G.grid(:, 1) > 0);
    I0 = (1-IF-IB);
    
    s = sF.*IF + sB.*IB;
    c = cF.*IF + cB.*IB + c0.*I0;
    u = param.u(c);

    % Construct finite-difference matrix for wealth evolution:
    Aa = FD_operator(G, s, zeros(G.J, 1), 1);

    % Solve linear system:
    A = Aa + Az;    
    B = (1/param.Delta + param.rho)*speye(G.J) - A;
    b = u + V / param.Delta;
    
    V_new = B\b;
    
    % Update value function:
    V_change = V_new - V;
    V = V_new;
    
    dist = max(max(abs(V_change)));
    if dist < param.crit, break; end
    if mod(iter, 1) == 0, fprintf('VFI: %.i    Remaining Gap: %.2d\n', iter, dist); end

end

end
