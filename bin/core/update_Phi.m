function phi = update_Phi(phi, grad_phi, alphai, options)
% UPDATE_PHI updates the yaw angle control inputs using gradient descent.

phi_lim = options.AMPC.Phi_lim;
Nc      = options.AMPC.Nc;
dJ      = grad_phi'; % Transpose gradient to match phi dimensions

% Simple gradient descent step
phis = phi(:, 1:Nc) - dJ * alphai;

% Apply hard constraints (clipping)
phis(phis > phi_lim(2)) = phi_lim(2);
phis(phis < phi_lim(1)) = phi_lim(1);

phi = [phis phis(:,end) * ones(1, length(phi) - Nc)];

end
