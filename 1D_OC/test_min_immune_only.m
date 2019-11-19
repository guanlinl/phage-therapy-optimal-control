% Testing the minimum level required to kill bacteria without phage
% within 72hrs.

Io = 8e6; % initial immune level
%% Fixed parameters 
% susceptible bacteria growth rate
p.r = 0.75;
% resistant bacteria growth rate
p.rp = 0.675;
% total bacteria carrying capacity
p.Kc = 1e10;
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

% parameters and initial conditions, wild-type = 1, Myd88 = 2;
parameter_choice = 2;
switch parameter_choice
    case 1 % wild-type
        p.Kc = 1e10; % total bacteria carrying capacity
        p.Ki = 2.4e7; % max capacity of immune response        
        % initial densities (fixed in this case)
        So = 7.4e7; Ro = 0; Po = 0; Io = 2.7e6; 
    case 2 % Myd88
        % total bacteria carrying capacity
        p.Kc = 8.5e11; 
        % initial densities
        So = 7.4e5; Ro = 0; Po = 0;
        % immue response is deactivated 
        p.Ki = Io;
end

%% No therapy sim
% time interval and fine grids parameters
% delayed therapy period
% run delay hours simulation to find the initial condition
delay_hours = 2;
p.t0 = 0; p.tf = delay_hours; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,1); % zero-injection rate

scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki]; 
ic_NoTherapy = [So;Ro;Po;Io]./scale_factor; 

p.P_in = 0;
x_state_NoTherapy = state_solver(ic_NoTherapy,p,u_NoTherapy); % simulation
t00 = p.t0:p.dt:p.tf;
initial_cond = x_state_NoTherapy(end, 2:end);

% rescale initial condition
y0 = initial_cond;
% redefine t
p.t0 = delay_hours; p.tf = 3*24; p.dt = 2e-4; p.N = round((p.tf - p.t0)/p.dt) + 1;
t = p.t0:p.dt:p.tf;
% rescale
q1 = 0;
t_start = 2; t_end = 3;
u_a = zeros(p.N,1); u_a(t >= t_start) = 1;
u_b = zeros(p.N,1); u_b(t >= t_end) = -1;
u_NC = q1.*(u_a + u_b); % naive control
x_state_NC = state_solver(y0,p,u_NC); % simulation

% stack data
x_state_NC = [x_state_NoTherapy; x_state_NC];
t = [t00, t];

% Population densities
figure(2); subplot(1,2,1); 
semilogy(t,x_state_NC(:,2).*scale_factor(1),'r','LineWidth',2); hold on;
semilogy(t,x_state_NC(:,3).*scale_factor(2),'b','LineWidth',2); hold on;
semilogy(t,x_state_NC(:,4).*scale_factor(3),'--r','LineWidth',2); hold on;
semilogy(t,x_state_NC(:,5).*scale_factor(4),'g','LineWidth',2); hold on;
legend('S, sensitive-bacteria','R, resistant-bacteria',...
        'P_{S}, phage','I, host immunity');
h = legend('S, sensitive-bacteria','R, resistant-bacteria',...
        'P_{S}, phage','I, host immunity');
set(h,'FontName','Times New Roman','FontSize',16);
legend boxoff;
axis([0,max(t),1e-2,1e10]); 
yticks([1e-2 1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',16);
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20); 
ylabel('Density', 'FontName', 'Times New Roman','FontSize',20); 
axis square;

% total bacteria population
figure(2); subplot(1,2,2); 
semilogy(t,(x_state_NC(:,2) + x_state_NC(:,3)).*...
    scale_factor(1),'k','LineWidth',3); 
ybars = [1e-2 1];
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)],...
    [0.6 0.6 0.6], 'FaceAlpha',.5);
legend('S + R, total bacteria');
h = legend('S + R, total bacteria');
set(h,'FontName','Times New Roman','FontSize',16);
legend boxoff;
axis([0,max(t),1e-2,1e10]); 
yticks([1e-2 1 1e2 1e4 1e6 1e8 1e10]); xticks(0: 12: p.tf); 
set(gca,'FontSize',16);
% ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xlabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20); 
ylabel('Density', 'FontName', 'Times New Roman','FontSize',20); 
axis square;


