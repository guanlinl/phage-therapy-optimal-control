
function result = treatment_pratical_multi_SDs(therapy_scenario,immune_level,u1, u2, ...
    t_seq, amp_fac)
%% ----------------------- Model parameters -----------------------
% Fixed parameters 
% susceptible bacteria growth rate
p.r = 0.75;
% resistant bacteria growth rate
p.rp = 0.675;
% total bacteria carrying capacity
p.Kc = 1e10;
% adsorption rate of phage:
p.phi = 5.4e-8;
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

p.P_in = 1e9; % ,maximal injection rate

% parameters and initial conditions, wild-type = 1, Myd88 = 2;
switch therapy_scenario
    case 1 % wild-type
        p.Kc = 1e10; % total bacteria carrying capacity
        p.Ki = 2.4e7; % max capacity of immune response
        % initial densities, R0 = 1
        So = 7.4e7; Ro = 1; Po = 0; Io = 2.7e6; P2o = 0; 
    case 2 % Myd88
        % total bacteria carrying capacity
        p.Kc = 8.5e11; 
        % initial densities, R0 = 1
        So = 7.4e5; Ro = 1; Po = 0; P2o = 0; Io = immune_level;
        % immue response is deactivated 
        p.Ki = Io;
end
p.dt = 2e-4;
scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki;p.Pc];  % scaled factor fixed

inj_data = pratical_injection_data(u1, u2, t_seq, p);
time_points = inj_data(1:5); total_dosage = inj_data(6:7).*(1 + amp_fac); % scaled, u

x_state = pratical_multi_SDs_solver(time_points, total_dosage, immune_level, p);

% compute total phage amount (dosage)
V_P = sum(total_dosage)*p.P_in; 

% total bacteria under phage control
N_sum = x_state(:,2).*scale_factor(1) + x_state(:,3).*scale_factor(2);
if isempty(find(N_sum < 1)) == 0
    result = [1, V_P]; % phage control works
else
    result = [0, 0]; % fail
end



end