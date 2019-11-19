% Numerical experiments: 2D_OC (parameter of interests: immune level)
% outputs: the population dynamics with some strategy
%          the plot of controls
%          the decay of cost function vs. iterations
%          the convergence indicator, \Theta vs. iteration

% parameter of interests
sim_case = 3;
switch sim_case
    case 1 % High immune basal level parameters
        Io = 8e6; % initial immune level
        weight_theta2 = 10;
        amp_fac = 0; % relaxed factor from 2D-OC to heuristic strategy
        Name1 = 'fig_2DOC_immune_high';
        Name2 = 'fig_2DOC_immune_high_practical';
    case 2 % Low immune basal level parameters
        Io = 3e6; % initial immune level
        weight_theta2 = 0.01; % read from text.file, not 0.03
        amp_fac = 2; % relaxed factor from 2D-OC to heuristic strategy, read from pratical_dosage_2D_SDs_immnue_level.txt
        Name1 = 'fig_2DOC_immune_low';
        Name2 = 'fig_2DOC_immune_low_practical';
    case 3 % very high immune 
        Io = 8.5e6; % initial immune level
        weight_theta2 = 1e11;
        amp_fac = 0; % relaxed factor from 2D-OC to heuristic strategy
        Name1 = 'fig_2DOC_immune_very_high';
        Name2 = 'fig_2DOC_immune_very_high_practical';
end 


%% Model parameters 
delay_hours = 2; % fixed
p.P_in = 1e9; % maximal injection rate of 2D control
% susceptible bacteria growth rate
p.r = 0.75;
% resistant bacteria growth rate
p.rp = 0.675;
% ------------------------------
% adsorption rate of phage, changes to make phage decay sensitive
% p.phi = 3.38e-8;
p.phi = 5.4e-8;
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
p.w = 0.07;
% maximum growth rate of immune response:
p.a = 0.97;
% conc. of bacteria at which imm resp growth rate is half its maximum:
p.Kn = 1e7;
% probability of emergence of phage-resistant mutation per cell division
p.m = 2.85e-8;

therapy_scenario = 2;
% parameters and initial conditions, wild-type = 1, Myd88 = 2;
switch therapy_scenario
    case 1 % wild-type
        p.Kc = 1e10; % total bacteria carrying capacity
        p.Ki = 2.4e7; % max capacity of immune response
        p.P_in = max_inj; % ,maximal injection rate
        % initial densities
        So = 7.4e7; Ro = 1; Po = 0; Io = 2.7e6; P2o = 0; 
    case 2 % Myd88
        % total bacteria carrying capacity
        p.Kc = 8.5e11; 
        % initial densities
        So = 7.4e5; Ro = 1; Po = 0; P2o = 0; 
        % immue response is deactivated 
        p.Ki = Io;
end

%% No therapy sim
% delayed therapy period
% run delay hours simulation to find the initial condition
p.t0 = 0; p.tf = delay_hours; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,2); % zero-injection rate

scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki;p.Pc]; 
ic_NoTherapy = [So;Ro;Po;Io;P2o]./scale_factor; 

x_state_NoTherapy = state_solver(ic_NoTherapy,p,u_NoTherapy); % simulation
initial_cond = x_state_NoTherapy(end, 2:end);
t00 = p.t0:p.dt:p.tf;


%% ----------------------- Control parameters -----------------------
% time interval and fine grids parameters
p.t0 = delay_hours; p.tf = 3*24; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
t = p.t0:p.dt:p.tf;
% Backtracking search parameters
p.beta_control = 0.5; p.alpha = 0.5; % step size parameters
p.theta0 = 1; p.theta1 = weight_theta2; p.theta2 = 1; % weights in cost functional

% initial guess of u(t)
initial_guess = 2;
switch initial_guess
    case 1 % random in [0,1]   
      u0 = rand(p.N,2)./2; % random initial guess
    case 2 
      u0 = 0.2.*ones(p.N,2); % uniformly = 0.5
    case 3 % step function
    t_prior = 10; 
    u0 = 1e-1.*ones(p.N,2); 
    % u0(find(t > t_prior),:) = 0.45;% no injection of phage
end
u = u0; 

% rescale initial states
y0 = initial_cond;

%% -------------------- Optimal control implementation ------------------
kk = 0;  N_iter = 30; tol_eps = 1e-5; % total iterations
% Hamiltonian-based Algorithm
for i = 1:N_iter
x_state = state_solver(y0,p,u);
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
% display the intermediate results
Performance = {'iterations';'(scaled) objective functional J';...
    '(scaled) convergence measure Theta'};
Numerics = [i; Ju; theta_u]; display(table(Performance,Numerics));

% break the OC loop
if abs(theta_u) < tol_eps
    break;
end
    
end

% compute total phage amount (dosage)
V_P_vec = sum(u).*p.P_in.*p.dt;
V_PS = V_P_vec(1); V_PR = V_P_vec(2); % two types of phage, total amount
V_P = V_PS + V_PR; % total total phage amount plus

% stack data
t = [t00,t];
x_state = [x_state_NoTherapy; x_state];
u = [u_NoTherapy;u];



%% Figures
% Population densities
figure(1); subplot(1,3,1); 
semilogy(t,x_state(:,2).*scale_factor(1),'r','LineWidth',2); hold on;
semilogy(t,x_state(:,3).*scale_factor(2),'b','LineWidth',2); hold on;
semilogy(t,x_state(:,4).*scale_factor(3),'--r','LineWidth',2); hold on;
semilogy(t,x_state(:,5).*scale_factor(4),'g','LineWidth',2); hold on;
semilogy(t,x_state(:,6).*scale_factor(5),'--b','LineWidth',2); hold on;
%legend('S, sensitive-bacteria','R, resistant-bacteria',...
%        'P_{S}, phage','I, host immunity','P_{R}, phage');
%h = legend('S, sensitive-bacteria','R, resistant-bacteria',...
%        'P_{S}, phage','I, host immunity','P_{R}, phage');
%set(h,'FontName','Times New Roman','FontSize',16);
%legend boxoff;
axis([0,max(t),1e-2,1e10]); 
yticks([1e-2 1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
ylabel('Density (/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
axis square;


figure(1); subplot(1,3,2); 
index_be = index_bacteria_elimination((x_state(:,2) + x_state(:,3)).*scale_factor(1),...
    1); % the index when the bacteria elimination
if isempty(index_be)
    semilogy(t,(x_state(:,2) + x_state(:,3)).*...
    scale_factor(1),'k','LineWidth',3); 
else
    semilogy(t(1:index_be),(x_state(1:index_be,2) + x_state(1:index_be,3)).*...
    scale_factor(1),'k','LineWidth',3); hold on;
    semilogy(t(index_be + 1:end),(x_state(index_be + 1:end,2) + x_state(index_be+1:end,3)).*...
        scale_factor(1),'--k','LineWidth',1.5); hold on;
    semilogy(t(index_be),(x_state(index_be,2) + x_state(index_be,3)).*scale_factor(1),'ok','MarkerSize',12);
    ybars = [1e-2 1];
    patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)],...
        [0.6 0.6 0.6], 'FaceAlpha',.5);
    %legend('S + R, total bacteria',...
    %    'S + R, total bacteria (extinct)');
    %h = legend('S + R, total bacteria',...
    %    'S + R, total bacteria (extinct)');
    %set(h,'FontName','Times New Roman','FontSize',16);
    %legend boxoff;
end
axis([0,max(t),1e-2,1e10]); 
yticks([1e-2 1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
ylabel('Density (/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
axis square;


% Optimal controls - injection strategies
figure(1); subplot(1,3,3); 
semilogy(t,u(:,1).*p.P_in + 1e-1,'r','LineWidth',2.5); hold on;
semilogy(t,u(:,2).*p.P_in + 1e-1,'b','LineWidth',2.5); hold on;
%semilogy(t,ones(p.N,1).*p.P_in,'--k','Linewidth',1);
axis([0,max(t),1,p.P_in*10]);
order_p = 0:2:log10(p.P_in*10);
yticks(10.^order_p); 
xticks(0: 12: p.tf); 
h_leg = legend('$\rho_{S}(t)$, phage $P_{S}$','$\rho_{R}(t)$, phage $P_{R}$','Location','south'); 
set(h_leg,'FontName', 'Times New Roman','FontSize',16, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
ylabel('Injection rate', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
axis square;
% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 40; Height = 10;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

print(ff, Name1, '-painters', '-depsc2','-r600');
print(ff, Name1, '-dpdf','-r600');


% -------------------------------------------------------------------------
%% Figures of multi-SDs
inj_data = pratical_injection_data(u(:,1), u(:,2), t, p);
time_points = inj_data(1:5); total_dosage = inj_data(6:7).*(1 + amp_fac);
x_state_SDs = pratical_multi_SDs_solver(time_points, total_dosage, Io, p);
x_state_SDs = pseudo_practical_process(x_state_SDs, p); % pseudo pratical conversion

% Population densities
figure(2); subplot(1,3,1); 
t = x_state_SDs(:,1);
semilogy(t,x_state_SDs(:,2).*scale_factor(1),'r','LineWidth',2); hold on;
semilogy(t,x_state_SDs(:,3).*scale_factor(2),'b','LineWidth',2); hold on;
semilogy(t,x_state_SDs(:,4).*scale_factor(3),'--r','LineWidth',2); hold on;
semilogy(t,x_state_SDs(:,5).*scale_factor(4),'g','LineWidth',2); hold on;
semilogy(t,x_state_SDs(:,6).*scale_factor(5),'--b','LineWidth',2); hold on;

legend('$S$, sensitive-bacteria','$R$, resistant-bacteria',...
        '$P_{S}$, phage','$I$, host immunity','$P_{R}$, phage', 'Position', [0.7 0.25 0.1 0.2]);
h = legend('$S$, sensitive-bacteria','$R$, resistant-bacteria',...
        '$P_{S}$, phage','$I$, host immunity','$P_{R}$, phage', 'Position',[0.7 0.25 0.1 0.2]);
set(h,'FontName', 'Times New Roman','FontSize',16, 'Interpreter','latex');
legend boxoff;

axis([0,max(t),1,1e10]); 
yticks([1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
ylabel('Density (/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
axis square;


figure(2); subplot(1,3,2); 
index_be = index_bacteria_elimination((x_state_SDs(:,2) + x_state_SDs(:,3)).*scale_factor(1),...
    1); % the index when the bacteria elimination
if isempty(index_be)
    semilogy(t,(x_state_SDs(:,2) + x_state_SDs(:,3)).*...
    scale_factor(1),'k','LineWidth',3); 
else
    semilogy(t(1:index_be),(x_state_SDs(1:index_be,2) + x_state_SDs(1:index_be,3)).*...
    scale_factor(1),'k','LineWidth',3); hold on;
    semilogy(t(index_be + 1:end),(x_state_SDs(index_be + 1:end,2) + x_state_SDs(index_be+1:end,3)).*...
        scale_factor(1),'--k','LineWidth',1.5); hold on;
    semilogy(t(index_be),(x_state_SDs(index_be,2) + x_state_SDs(index_be,3)).*scale_factor(1),'ok','MarkerSize',12);
    ybars = [1e-2 1];
    patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)],...
        [0.6 0.6 0.6], 'FaceAlpha',.5);
    legend('$S + R$, total bacteria',...
        '$S + R$, total bacteria (extinct)', 'Position', [0.73 0.55 0.1 0.1]);
    h = legend('$S + R$, total bacteria',...
        '$S + R$, total bacteria (extinct)', 'Position', [0.73 0.55 0.1 0.1]);
    set(h,'FontName', 'Times New Roman','FontSize',16, 'Interpreter','latex');
    legend boxoff;
end
axis([0,max(t),1,1e10]); 
yticks([1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
ylabel('Density (/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex'); 
axis square;


% output figures in eps 
fff = figure(2);
fff.Units = 'inches';
Width = 40; Height = 10;

fff.PaperSize = [Width, Height];
fff.PaperPosition = [0 0 Width, Height];
fff.Position = [0 0 Width, Height];

print(fff, Name2, '-painters', '-depsc2','-r600');
print(fff, Name2, '-dpdf','-r600');




