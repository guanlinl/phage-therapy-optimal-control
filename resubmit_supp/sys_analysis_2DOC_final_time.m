%% Explore the 2D-OC strategy vs. final time variation
%  immune deficient case

% define final time set from 2 days to 5 days
tf_set = 48:3:96;
immune_level = 6e6; % fix intermediate level 

% fixed parameters
therapy_scenario = 2; %  immune efficient case
phage_absorp = 5.4e-8; % fixed absorption rate, high
phage_decay = 0.07; % fixed

% adaptive search of theta2, initial grid
power_theta_range_init = [-11, 0, 11]; 
num_adpt = 10; % number of search iterations

% initialization of data saver
minimal_dosage_data_2D_OC = zeros(2, length(tf_set)); % for each time, we need two phage compositions
theta2_2D_OC = zeros(1, length(tf_set)); % record tunned parameters
inj_time_data_2D_OC = zeros(2, length(tf_set));

% 2D-OC dosage data by varying theta2
for k = 1:length(tf_set)
    FinalTime = tf_set(k);
   
    % sweeping on initial mesh
    result_a = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(1), FinalTime);
    result_b = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(2), FinalTime);
    result_c = treatment_2D_OC(therapy_scenario, immune_level,...
                  10.^power_theta_range_init(3), FinalTime);
    
              
    if result_a(1) <= 0 % nothing works
        minimal_dosage_data_2D_OC(k) = 0;
        theta2_2D_OC(k) = 0; % default ...
    elseif result_c(1) > 0 % everything works
        minimal_dosage_data_2D_OC(1, k) = result_c(3);
        minimal_dosage_data_2D_OC(2, k) = result_c(4);
        inj_time_data_2D_OC(1, k) = result_c(5);
        inj_time_data_2D_OC(2, k) = result_c(6);
        
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
                                10.^power_theta_curr, FinalTime);     
            if result_curr(1) > 0 % works
                power_theta_cand_left = power_theta_curr;
            else % not work
                power_theta_cand_right = power_theta_curr;
            end
            l = l + 1; % update iteration num
        end
        minimal_dosage_data_2D_OC(1, k) = result_curr(3);
        minimal_dosage_data_2D_OC(2, k) = result_curr(4);
        inj_time_data_2D_OC(1, k) = result_curr(5);
        inj_time_data_2D_OC(2, k) = result_curr(6);
        theta2_2D_OC(k) = 10^power_theta_curr;
    end
    
    disp(k);
end

% write data file
y = [tf_set; minimal_dosage_data_2D_OC; inj_time_data_2D_OC; theta2_2D_OC];
fileID = fopen('minimal_dosage_data_2D_OC_tf.txt','w');
fprintf(fileID,...
    'tf_set, minimal_dosage_PS, minimal_dosage_PR, ts, tr, theta2 \n\n');
fprintf(fileID,'%f %f %f %f %f %f\n',y);
fclose(fileID);



