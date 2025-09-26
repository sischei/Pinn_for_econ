function hjb = savings(V, G, p)

num0 = 1e-8; % numerical 0 for upwind scheme

VaF = zeros(G.J, p.Nd);
VaB = zeros(G.J, p.Nd);
for j = 1:p.Nd
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

VaF(G.grid(:, 1) == 1, :) = p.u1(G.income(G.grid(:, 1) == 1, :));
VaB(G.grid(:, 1) == 0, :) = p.u1(G.income(G.grid(:, 1) == 0, :));

% Consumption:
cF = p.u1inv(VaF);
cB = p.u1inv(VaB);
c0 = G.income;

% Savings:
sF = G.income - cF;
sB = G.income - cB;

% Upwinding:
IF = (sF > num0);        % BC takes care of this: (G.grid(:, 1) < 1)
IB = (sB <-num0) & ~IF;  % BC takes care of this: (G.grid(:, 1) > 0)
I0 = ~IF & ~IB;

hjb.s = sF.*IF + sB.*IB;
hjb.c = cF.*IF + cB.*IB + c0.*I0;
hjb.u = p.u(hjb.c);

end

