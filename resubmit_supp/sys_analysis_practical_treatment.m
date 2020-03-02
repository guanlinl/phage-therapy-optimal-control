
treatment_data_2D_OC = importdata('minimal_dosage_data_2D_OC_tf.txt');
tf_set = treatment_data_2D_OC.data(:,1);

% minimal dosage data read
dosage_PS = treatment_data_2D_OC.data(:,2);
% minimal dosage data read
dosage_PR = treatment_data_2D_OC.data(:,3);

%
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
therapy_scenario = 2;
switch therapy_scenario
    case 1 % wild-type
        p.Kc = 1e10; % total bacteria carrying capacity
        p.Ki = 2.4e7; % max capacity of immune response
        % initial densities
        So = 7.4e7; Ro = 1; Po = 0; Io = 2.7e6; P2o = 0; 
    case 2 % Myd88
        % total bacteria carrying capacity
        p.Kc = 8.5e11; 
        % initial densities
        So = 7.4e5; Ro = 1; Po = 0; P2o = 0; Io = 6e6;
        % immue response is deactivated 
        p.Ki = Io;
end

dtt = 1e-3;

% delayed therapy period
delay_hours = 2;
% run delay hours simulation to find the initial condition
p.t0 = 0; p.tf = delay_hours; p.dt = dtt; p.N = round((p.tf - p.t0)/p.dt) + 1;
u_NoTherapy = zeros(p.N,2); % zero-injection rate

scale_factor = [p.Kd;p.Kd;p.Pc;p.Ki;p.Pc]; 
ic_NoTherapy = [So;Ro;Po;Io;P2o]./scale_factor; 

x_state_NoTherapy = state_solver(ic_NoTherapy,p,u_NoTherapy); % simulation
initial_cond = x_state_NoTherapy(end, 2:end);

dosage_PS_prac = zeros(length(tf_set), 1); dosage_PR_prac = zeros(length(tf_set), 1); 

p.increment = 1e-2;
for t_id = 1:length(tf_set)
    % start treatment 
    p.t0 = delay_hours; p.tf = tf_set(t_id); % for numerical purpose
    p.dt = dtt; p.N = round((p.tf - p.t0)/p.dt) + 1;
    t = p.t0:p.dt:p.tf;
    u_NoTherapy = zeros(p.N,2); % zero-injection rate
    
    eliminate_flag = false; kk = 0;
    
    while eliminate_flag == false 
        ic_therapy = initial_cond + ...
            [0, 0, dosage_PS(t_id)/p.Pc, 0, dosage_PR(t_id)/p.Pc].*(1 + kk*p.increment);
        x_state_therapy = state_solver(ic_therapy,p,u_NoTherapy); % simulation

        % total bacteria under phage control
        N_sum = x_state_therapy(:,2).*scale_factor(1) + x_state_therapy(:,3).*scale_factor(2);
        if isempty(find(N_sum < 1)) == 0
            eliminate_flag = true;
        else
            kk = kk + 1; 
        end
    end
    display(kk)
    dosage_PS_prac(t_id) = dosage_PS(t_id)*(1 + kk*p.increment);
    dosage_PR_prac(t_id) = dosage_PR(t_id)*(1 + kk*p.increment);
end


% write data file
y = [tf_set'; [dosage_PS_prac'; dosage_PR_prac']];
fileID = fopen('minimal_dosage_data_2D_OC_tf_practical.txt','w');
fprintf(fileID,...
    'tf_set, minimal_dosage_PS, minimal_dosage_PR \n\n');
fprintf(fileID,'%f %f %f\n',y);
fclose(fileID);





