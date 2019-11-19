% Results of 1D-OC phage therapy
% Input:
% therapy_scenario - immune efficient (1) or immune deficient (2)
% immune_level - initial immune response level
% phage_absorp - phage adsorption rate
% phage_decay - phage decay rate
% weight_theta2 - weight of theta2, i.e., the regularizartion of phage
%                 injection rate
% Output: 
% [0, 0] - if control doesn't work
% [1, total_dosage] - if control works (bacterial elimination)

function result = treatment_1D_OC(therapy_scenario, immune_level,...
    phage_absorp, phage_decay,  weight_theta2)

%% ----------------------- Model parameters -----------------------
% Fixed parameters 
% susceptible bacteria growth rate
p.r = 0.75;
% resistant bacteria growth rate
p.rp = 0.675;
% total bacteria carrying capacity
p.Kc = 1e10;
% ------------------------------
% adsorption rate of phage, changes to make phage decay sensitive
p.phi = phage_absorp;
% ------------------------------
% phage density at half saturation
p.Pc = 1.5e7;
% immune response killing rate parameter:
p.ep = 8.2e-8;
% bacterial conc. at which immune response is half as effective:
p.Kd = 4.1e7;
% burst size of phage:
p.beta = 100;
% decay rate of phage:
p.w = phage_decay;
% maximum growth rate of immune response:
p.a = 0.97;
% conc. of bacteria at which imm resp growth rate is half its maximum:
p.Kn = 1e7;
% probability of emergence of phage-resistant mutation per cell division
p.m = 2.85e-8;

% parameters and initial conditions, wild-type = 1, Myd88 = 2;
switch therapy_scenario
    case 1 % wild-type
        p.Kc = 1e10; % total bacteria carrying capacity
        p.Ki = 2.4e7; % max capacity of immune response
        % initial densities (fixed in this case)
        So = 7.4e7; Ro = 1; Po = 0; Io = 2.7e6;
    case 2 % Myd88
        % total bacteria carrying capacity
        p.Kc = 8.5e11; 
        % initial densities
        So = 7.4e5; Ro = 1; Po = 0; 
        Io = immune_level;
        % immue response is deactivated 
        p.Ki = Io;
end

p.P_in = 1e9; % ,maximal injection rate

% run delay hours simulation to find the initial condition
delay_hours = 2; % delayed therapy period, 2hrs
p.t0 = 0; p.tf = delay_hours; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,1); % zero-injection rate

scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki]; 
ic_NoTherapy = [So;Ro;Po;Io]./scale_factor; 

x_state_NoTherapy = state_solver(ic_NoTherapy,p,u_NoTherapy); % simulation
initial_cond = x_state_NoTherapy(end, 2:end);

% ------------------------------------------------------------------
%% ----------------------- Control parameters -----------------------
% time interval and fine grids parameters
p.t0 = delay_hours; p.tf = 24*3; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
t = p.t0:p.dt:p.tf;
% Backtracking search parameters
p.beta_control = 0.5; p.alpha = 0.5; % step size parameters

p.theta0 = 1; p.theta1 = weight_theta2; p.theta2 = 1; % weights in cost functional


% initial guess of u(t), 1D
initial_guess = 2;
switch initial_guess
    case 1 % random in [0,1]
      u0 = rand(p.N,1); % random initial guess
    case 2 
      u0 = 0.5.*ones(p.N,1); % uniformly = 0.5
    case 3 % step function
    u0 = 1e-1.*ones(p.N,1); 
    % t_prior = 10; u0(find(t > t_prior)) = 1;% no injection of phage
end
u = u0; 

% rescale initial condition
y0 = initial_cond;

%% -------------------- Optimal control implementation ------------------
N_iter = 20; tol_eps = 1e-10; % total iterations and convergence tolerance

for i = 1:N_iter
x_state = state_solver(y0,p,u); % forward ode solver 
Ju = cost_obj(x_state, u, p); % compute cost J[u]
co_state = costate_solver(x_state,p); % compute co-state
% optimal policy
u_star = optimal_u(co_state, p);
% Hamiltonian
H_u_star = Hamiltonian(x_state, co_state, u_star, p);
H_u_temp = Hamiltonian(x_state, co_state, u, p);
% Theta_measure of u;
theta_u = sum(H_u_star - H_u_temp)*p.dt;
% compute kk = 0 case
x_state_temp = state_solver(y0,p,u_star); 
J_temp = cost_obj(x_state_temp, u_star, p);
Delta_J = J_temp - Ju; kk = 0;
 
    while (Delta_J >= p.alpha*p.beta_control^kk*theta_u) % step size check
        kk = kk + 1;
        u_curr = (1 - p.beta_control^kk).*u +  (p.beta_control^kk).*u_star;
        x_state_temp = state_solver(y0,p,u_curr);
        J_temp = cost_obj(x_state_temp, u_curr, p);
        Delta_J = J_temp - Ju;
        if abs(Delta_J) < 1e-13 
            break; 
        end
    end
u = (1 - p.beta_control^kk).*u + (p.beta_control^kk).*u_star; % update control 

% break the OC loop
if abs(theta_u) < tol_eps
    break;
end

end

% compute total phage amount (dosage)
V_P = sum(u).*p.P_in.*p.dt; % total total phage amount plus

% total bacteria under phage optimal control
N_sum = x_state(:,2).*scale_factor(1) + x_state(:,3).*scale_factor(2);
if isempty(find(N_sum < 1)) == 0
    result = [1, V_P]; % phage OC works
else
    result = [0, 0]; % fail
end

end