% ODE solution of multi-single dose 
% INPUT: 
% time points: [0, no_therapy_end, tS_inj, tR_inj, therapy_end]
% injection amount: total dosage of u_S, total dosage of u_R, scaled
% parameter information: immune density, others are fixed.

% Output: solutions of ODEs, [time, state_solutions]

function results = pratical_multi_SDs_solver(time_points, injection_amount, immune_level, p)

% assign parameter of interests
Io = immune_level; % initial immune level
P_S_scaled_ic = injection_amount(1)*p.P_in/p.Pc; % scaled by p.Pc, 1 is the sum of u input
P_R_scaled_ic = injection_amount(2)*p.P_in/p.Pc; % scaled by p.Pc

t_init = time_points(1); % default: t_init  = 0
delay_hours = time_points(2);% default: delay_hours  = 2
tS_inject = time_points(3); 
tR_inject = time_points(4);
t_end = time_points(5); % default: t_end  = 72

if tS_inject == tR_inject % boundary case, not likely to happen
    tR_inject = tR_inject + p.dt; % arbitrarily seperate injection timings
end
    
% initial condition added
if tS_inject < tR_inject % inject PS first then PR
    t1 = tS_inject; t2 = tR_inject;
    ic_t1_added = [0, 0, P_S_scaled_ic, 0, 0];
    ic_t2_added = [0, 0, 0, 0, P_R_scaled_ic];
else % inject PR first then PS
    t1 = tR_inject; t2 = tS_inject;
    ic_t1_added = [0, 0, 0, 0, P_R_scaled_ic];
    ic_t2_added = [0, 0, P_S_scaled_ic, 0, 0];
end

    
% ---------------------- Model parameters (fixed!!!) ------------------
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
% -----------------------------------------------------------------------

p.dt = 2e-4; % time step is fixed
scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki;p.Pc];  % scaled factor fixed

% ----------------------- no therapy simulation ------------------------
% define initial and end times in this case, t1 >= delayed_hours
p.t0 = t_init; p.tf = t1;  p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,2); % zero-injection rate

ic_NoTherapy = [So;Ro;Po;Io;P2o]./scale_factor; % ic at initial, t = 0

% no therapy, x_state_NoTherapy = [time, x_sol]
x_state_NoTherapy = state_solver(ic_NoTherapy,p,u_NoTherapy); % simulation

% final state, from case 0 to case 1
ic_01 = x_state_NoTherapy(end,2:end);
ic1 = ic_01 + ic_t1_added;
% -----------------------------------------------------------------------


%% ----------------------- after first phage shot ------------------------
% define initial and end times in this case
p.t0 = t1; p.tf = t2;  p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,2); % zero-injection rate
x_state_FirstShot = state_solver(ic1,p,u_NoTherapy); % simulation

% final state, from case 1 to case 2
ic_12 = x_state_FirstShot(end,2:end);
ic2 = ic_12 + ic_t2_added;
% -----------------------------------------------------------------------


%% ----------------------- after second phage shot ------------------------
% define initial and end times in this case
p.t0 = t2; p.tf = t_end;  p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,2); % zero-injection rate
x_state_SecondShot = state_solver(ic2,p,u_NoTherapy); % simulation
% -------------------------------------------------------------------------

results = [x_state_NoTherapy; x_state_FirstShot; x_state_SecondShot]; 

end








