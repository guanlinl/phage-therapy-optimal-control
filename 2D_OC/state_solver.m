% Numerical solvers for ODEs systems
% Input : state equations (and fixed parameters),
%         initial conditions, time interval, the external control signals u(t)
%         method, explicit Euler (1), RK4 (2); maybe more later (implicit schemes ?)
% Output: state solutions with uniformed time grids

function state_solution = state_solver(ic,p,u_input)
dt = p.dt; N = p.N; % number of grids
t = p.t0:dt:p.tf; % time grids [t0, t0 + dt, ..., t0 + N*dt];
solution_saver = zeros(5,N); % initialize saver
solution_saver(:,1) = ic;  % initial condition, ic is a column vec

%% -------------------- Numerical solvers ---------------------
% forward euler method (explicit solver)
        for i = 1:N-1
            f = rpsODE(t(i),solution_saver(:,i),p,u_input(i,:)); % dxdt
            solution_saver(:,i+1) = solution_saver(:,i) + f.*dt;
        end
state_solution = [t' solution_saver'];

end
