function output_time_series = pseudo_practical_process(pop_data, p)
output_time_series = pop_data(:,1); % stack time index first
for k = 2:length(pop_data(1,:))
    if k == 2 || k == 3 % two bacteria populations
        temp_ind = min(find(pop_data(:,k) < 1/p.Kd));
        pop_data(temp_ind:end,k) = 0;        
    end
    output_time_series = [output_time_series, pop_data(:,k)];
end
end