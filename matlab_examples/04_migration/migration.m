function [pi, Api, V_gain] = migration(V, p)

Vmig = reshape(V, [p.Naz, p.Nr]);
E_exp_Vmig = (p.etk * exp(p.theta * Vmig'))';

% Migration shares:
exp_Vmig = exp(p.theta * Vmig);
numerator = permute(repmat(exp_Vmig, [1, 1, p.Nr]), [1 3 2]) .* ...
            permute(repmat(p.etk, [1 1 p.Naz]), [3 1 2]);
pi = reshape(numerator ./ repmat(E_exp_Vmig, [1 1 p.Nr]), [p.Ntot p.Nr]);

% Indices for migration matrix:
x = repmat( (1:p.Ntot)' , [p.Nr,1] );
y = repmat( (1:p.Naz)'  , [p.Nr,1] ) + p.Naz * (0:(p.Nr-1)); 
y = y(:);

% Convert index pairs to linear indices:
mig_ind = sub2ind([p.Ntot, p.Ntot], x, y);

% Migration matrix:
mig = sparse(p.Ntot, p.Ntot);
mig(mig_ind) = pi(:);
Api = p.mu * (mig - speye(p.Ntot));

% Selection welfare gain:
kappa_bar = kron(p.kappa, ones(p.Naz, 1)) - 1/p.theta*log(pi);
V_gain = sum(pi.*kappa_bar, 2);

end

