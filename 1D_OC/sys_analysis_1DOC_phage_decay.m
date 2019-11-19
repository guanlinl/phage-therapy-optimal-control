%% Explore the 1D-OC strategy and SD strategy vs. phage decay rate
%  immune efficient case

% define ranges and scenario
exp_w = -2:0.2:1; w_set = 10.^exp_w; % phage decay range, 1e-2 ~ 1e2

therapy_scenario = 1; %  immune efficient case
phage_absorp = 3.38e-8; % fixed absorption rate, low
immune_level = 0; % in the case of immune efficient, default

% adaptive search of theta2, initial grid
power_theta_range_init = [-11, 0, 11]; 
num_adpt = 8; % number of search iterations

% initialization of data saver
minimal_dosage_data_1D_OC = zeros(1, length(w_set));
theta2_1D_OC = zeros(1, length(w_set));

% 1D-OC dosage data by varying theta2
for k = 1:length(w_set)
    phage_decay = w_set(k);
    
    % sweeping on initial mesh
    result_a = treatment_1D_OC(therapy_scenario, immune_level,...
                  phage_absorp, phage_decay,  10.^power_theta_range_init(1));
    result_b = treatment_1D_OC(therapy_scenario, immune_level,...
                  phage_absorp, phage_decay,  10.^power_theta_range_init(2));
    result_c = treatment_1D_OC(therapy_scenario, immune_level,...
                  phage_absorp, phage_decay,  10.^power_theta_range_init(3));
    
              
    if result_a(1) <= 0 % nothing works
        minimal_dosage_data_1D_OC(k) = 0;
        theta2_1D_OC(k) = 0; % default ...
    elseif result_c(1) > 0 % everything works
        minimal_dosage_data_1D_OC(k) = result_c(2);
        theta2_1D_OC(k) = 1e11;        
    else % result_a(1) > 0, result_c(2) < 0, middle
        % define left and right values of bi-search
        if result_b(1) > 0 % only c doesn't work 
           power_theta_cand_left = power_theta_range_init(2);
           power_theta_cand_right = power_theta_range_init(3);
        else % only a works
           power_theta_cand_left = power_theta_range_init(1);
           power_theta_cand_right = power_theta_range_init(2);
        end
        l = 0; result_curr = [0,0]; % initialize 
        while l < num_adpt || result_curr(1) <= 0 % break if l >= num_adp and result_curr(1) > 0
            power_theta_curr = (power_theta_cand_left + ...
                               power_theta_cand_right)/2;
            result_curr = treatment_1D_OC(therapy_scenario, immune_level,...
                  phage_absorp, phage_decay,  10.^power_theta_curr);     
            if result_curr(1) > 0 % works
                power_theta_cand_left = power_theta_curr;
            else % not work
                power_theta_cand_right = power_theta_curr;
            end
            l = l + 1; % update iteration num
        end
        minimal_dosage_data_1D_OC(k) = result_curr(2);
        theta2_1D_OC(k) = 10^power_theta_curr;
    end
    
    % check if needs to continue to the next phage decay rate
    if minimal_dosage_data_1D_OC(k) == 0 % doesn't exist curative treatment
        break;
    end
    disp(k);
end


% write data file
y = [w_set; minimal_dosage_data_1D_OC; theta2_1D_OC];
fileID = fopen('minimal_dosage_data_1D_OC_phage_decay.txt','w');
% fprintf(fileID, 'minimal_dosage_data_1D_OC vs. phage_decay rate \n\n');
fprintf(fileID,'%f %f %f\n',y);
fclose(fileID);



