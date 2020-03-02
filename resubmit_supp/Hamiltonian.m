function H = Hamiltonian(x_sol, lambda_sol, u_sol, p)
% H(x, u, p) = <p, f(x,u)> + L(x,u)
% ---------- balancing weights -------------
theta0 = p.theta0; theta1 = p.theta1; theta2 = p.theta2;
% ------ state, costate and controls -------
x1 = x_sol(:,2); x2 = x_sol(:,3); x3 = x_sol(:,4); 
x4 = x_sol(:,5); x5 = x_sol(:,6);
lam1 = lambda_sol(:,2); lam2 = lambda_sol(:,3); 
lam3 = lambda_sol(:,4); lam4 = lambda_sol(:,5);
lam5 = lambda_sol(:,6);
u1 = u_sol(:,1); u2 = u_sol(:,2);

%----------- Loading parameters-----------
% susceptible bacteria growth rate
r = p.r;
% resistant bacteria growth rate
rp = p.rp;
% total bacteria carrying capacity
Kc = p.Kc;
% adsorption rate of phage:
phi = p.phi;
% phage density at half saturation
Pc = p.Pc;
% immune response killing rate parameter:
ep = p.ep;
% bacterial conc. at which immune response is half as effective:
Kd = p.Kd;
% burst size of phage:
beta = p.beta;
% decay rate of phage:
w = p.w;
% maximum growth rate of immune response:
a = p.a;
% max capacity of immune response:
Ki = p.Ki;
% conc. of bacteria at which imm resp growth rate is half its maximum:
Kn = p.Kn;
% probability of emergence of phage-resistant mutation per cell division
m = p.m;
% phage injection amount is accomplished in 1hr
P_in = p.P_in; 
%----------- rescaled parameters-----------
% recaled parameters will be noted by '(.)0' notation
ep0 = ep*Ki; 
q0 = P_in/Pc; 
Kcd0 = Kc/Kd; 
Knd0 = Kn/Kd; 
Kpd0 = Pc/Kd;
phi0 = phi*Kd;
% infection rate function (PS-R model)
F_P1 = phi0.*(x3./(1 + x3));
F_P2 = phi0.*(x5./(1 + x5));

n = x1 + x2;  % total (recaled) bacteria

% RHS of systems of equations
f1 = r.*x1.*(1 - n./Kcd0).*(1-m) - Kpd0.*x1.*F_P1 - ep0.*x1.*x4./(1 + n);
f2 = rp.*x2.*(1 - n./Kcd0) + r.*x1.*(1 - n./Kcd0).*m - ep0.*x2.*x4./(1 + n) - ...
    Kpd0.*x2.*F_P2;
f3 = beta.*x1.*F_P1 - phi0.*x1.*x3 - w.*x3 + q0.*u1;
% Change immune response
f4 = a.*x4.*(n./(n + Knd0)).*(1 - x4);
% Change in scaled phage 2 population
f5 = beta.*x2.*F_P2 - phi0.*x2.*x5 - w.*x5 + q0.*u2;
    
H = (f1.*lam1 + f2.*lam2 + f3.*lam3 + f4.*lam4 + f5.*lam5 + ...
    theta0.*(x1 + x2))' + (theta1/2).*vecnorm(u_sol',2).^2;
end


