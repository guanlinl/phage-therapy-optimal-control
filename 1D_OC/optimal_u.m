function u_star = optimal_u(lambda_sol,p)
q0 = p.P_in/p.Pc; 
% dependent costate
lam3 = lambda_sol(:,4); % lambda_3

% ----------- maximum principle -----------
u_star = -q0.*lam3./p.theta1; 
u_star(find(lam3 < -p.theta1/q0)) = 1; 
u_star(find(lam3 >= 0)) = 0;
end