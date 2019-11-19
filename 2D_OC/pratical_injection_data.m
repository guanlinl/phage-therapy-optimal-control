% heuristic approach of multi-SDs given 2D-OC sols
% INPUT: injection strategies given by 2D-OCs
% OUTPUT: [0, no_therapy_end, tS_inj, tR_inj, therapy_end, total_uS, total_uR]

function results = pratical_injection_data(uS, uR, t_seq, p)
max_uS_seq = find(uS == max(uS)); time_max_uS = t_seq(max_uS_seq);
max_uR_seq = find(uR == max(uR)); time_max_uR = t_seq(max_uR_seq);
tS_inj = time_max_uS(1); tR_inj = time_max_uR(1);

total_uS = sum(uS)*p.dt; total_uR = sum(uR)*p.dt;

results = [0, 2, tS_inj, tR_inj, 72, total_uS, total_uR];

end

