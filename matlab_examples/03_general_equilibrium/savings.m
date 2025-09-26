function policies = savings(V, G, p)

num0 = 1e-8; % numerical 0 for upwind scheme

% Forward and backward approximations:
VaF = zeros(G.J, p.discrete_types);
VaB = zeros(G.J, p.discrete_types);
for j = 1:p.discrete_types
    VaF(:, j) = deriv_sparse(G, V(:, j), 1, 'D1F', num2str(j));
    VaB(:, j) = deriv_sparse(G, V(:, j), 1, 'D1B', num2str(j));
end

if max(max(abs( VaF(G.grid(:, 1) == 1, :) - ...
    p.u1(G.income(G.grid(:, 1) == 1, :)) ))) > 1e-8
    fprintf('Right boundary off \n');
end
if max(max(abs( VaB(G.grid(:, 1) == 0, :) - ...
    p.u1(G.income(G.grid(:, 1) == 0, :)) ))) > 1e-8
    fprintf('Left boundary off \n');
end

% VaF(G.grid(:, 1) == 1, :) = p.u1(G.income(G.grid(:, 1) == 1, :));
% VaB(G.grid(:, 1) == 0, :) = p.u1(G.income(G.grid(:, 1) == 0, :));

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

policies.s = sF.*IF + sB.*IB;
policies.c = cF.*IF + cB.*IB + c0.*I0;
policies.u = p.u(policies.c);


% Tests:
assert(all(all( policies.s(G.grid(:, 1) == 1,:) <= 0 )));
assert(all(all( policies.s(G.grid(:, 1) == 0,:) >= 0 )));

end

