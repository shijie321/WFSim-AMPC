function Phi = update_Phi(phi, grad_phi, alphai, options)
% UPDATE_PHI updates the yaw angle control inputs using gradient descent.

phi_lim = options.AMPC.Phi_lim;
Nc      = options.AMPC.Nc;
dJ      = grad_phi'; % Transpose gradient to match phi dimensions

% Simple gradient descent step
phi = phi - dJ * alphai;

% Apply hard constraints (clipping)
phi(phi > phi_lim(2)) = phi_lim(2);
phi(phi < phi_lim(1)) = phi_lim(1);

end