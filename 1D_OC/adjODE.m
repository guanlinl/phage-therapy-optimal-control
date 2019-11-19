% adjoint equations (costate eqs)
% Input :  state eqs (forward step state) : x(t_i+1)
%          adjoint eqs (forward step state) : lam(t_i+1)
%          control variable u(t_i+1), (may not be used)
%          time t_i+1, (may not be used)
% Output : the differential vector

function dlamdt = adjODE(t,x,lam,p)
% ---------- state variables ---------------
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); 
% ---------- Grad_x[L(x,u)] -----------------
grad_L = [p.theta0; p.theta0; 0; 0];

% ---------- Jacobian matrix of x' = f(x,u) --------
J = Jacobian_mat(x1,x2,x3,x4,p); 

% ---------- construct the costate diff -----------
dlamdt = - (J')*lam - grad_L; % return the column vec
end

