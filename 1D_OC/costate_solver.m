% Numerical solvers for costate
% Input : state equations (and fixed parameters),
%         terminal conditions, time interval, the external control signals u(t)
%         method, backward Euler (1), RK4 (2); maybe more later (implicit schemes ?)
% Output: costate solutions with uniformed time grids

function costate_solution = costate_solver(state_sol,p)
dt = p.dt; N = p.N; % number of grids
t = p.t0:dt:p.tf; % time grids [t0, t0 + dt, ..., t0 + N*dt];

lam_sol = zeros(4,N); % initialization
lam_terminal = [p.theta2;p.theta2;0;0]; 
lam_sol(:,N) = lam_terminal; % terminal condition

state_sol = state_sol(:,2:5); state_sol = state_sol';
%% -------------------- Numerical solvers ---------------------
 % backward euler method 
   for i = 1:N-1
        f_adj = adjODE(t(N + 1 - i),state_sol(:,N + 1 - i),lam_sol(:,N + 1 - i),p); % dlambda/dt
        lam_sol(:,N-i) = lam_sol(:, N - i + 1) - f_adj.*dt; % backward euler
   end
            
costate_solution = [t' lam_sol'];

end
