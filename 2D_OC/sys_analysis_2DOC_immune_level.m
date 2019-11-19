%% Explore the 2D-OC strategy vs. immune level
%  immune deficient case

% define ranges and scenario
immune_level_set = 3e6:5e5:8.5e6; % immune level 3e6 ~ 8e6

therapy_scenario = 2; %  immune efficient case
phage_absorp = 5.4e-8; % fixed absorption rate, high
phage_decay = 0.07; % fixed

% adaptive search of theta2, initial grid
power_theta_range_init = [-11, 0, 11]; 
num_adpt = 8; % number of search iterations

% initialization of data saver
minimal_dosage_data_2D_OC = zeros(1, length(immune_level_set));
theta2_2D_OC = zeros(1, length(immune_level_set));

% 1D-OC dosage data by varying theta2
for k = 1:length(immune_level_set)
    immune_level = immune_level_set(k);
          
    % sweeping on initial mesh
    result_a = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(1));
    result_b = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(2));
    result_c = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(3));
    
              
    if result_a(1) <= 0 % nothing works
        minimal_dosage_data_2D_OC(k) = 0;
        theta2_2D_OC(k) = 0; % default ...
    elseif result_c(1) > 0 % everything works
        minimal_dosage_data_2D_OC(k) = result_c(2);
        theta2_2D_OC(k) = 1e11;
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
            result_curr = treatment_2D_OC(therapy_scenario, immune_level,...
                                10.^power_theta_curr);     
            if result_curr(1) > 0 % works
                power_theta_cand_left = power_theta_curr;
            else % not work
                power_theta_cand_right = power_theta_curr;
            end
            l = l + 1; % update iteration num
        end
        minimal_dosage_data_2D_OC(k) = result_curr(2);
        theta2_2D_OC(k) = 10^power_theta_curr;
    end
    
    
    disp(k);
end

% write data file
y = [immune_level_set; minimal_dosage_data_2D_OC; theta2_2D_OC];
fileID = fopen('minimal_dosage_data_2D_OC_immnue_level.txt','w');
fprintf(fileID, 'immnue level, minimal_dosage_data_2D_OC, theta2 \n\n');
fprintf(fileID,'%f %f %f\n',y);
fclose(fileID);



