function u_star = optimal_u(lambda_sol,p)
q0 = p.P_in/p.Pc; 
% dependent costate
lam3 = lambda_sol(:,4); % lambda_3
lam5 = lambda_sol(:,6); % lambda_5
lam_vec = [lam3, lam5];
% ----------- maximum principle -----------
u_hat = -(q0/p.theta1).*lam_vec; 
u_star = zeros(size(u_hat));
for k = 1:p.N 
    u_star(k,:) = proj_constrain_set(u_hat(k,:)')';
end
end