function cc = cost_obj(x_sol, u_sol, p)
% ---------- balancing weights -------------
theta0 = p.theta0; theta1 = p.theta1; theta2 = p.theta2;
x1 = x_sol(:,2); x2 = x_sol(:,3); x3 = x_sol(:,4); x4 = x_sol(:,5); 
% ----------- cost ------------------------
cc = sum(theta0.*(x1 + x2) + ((theta1)/2).*(u_sol.^2))*p.dt + ...
    theta2*(x1(end) + x2(end));
end



