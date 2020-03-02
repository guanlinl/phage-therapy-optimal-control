function index_be = index_bacteria_elimination(density, threhold)
index_set = find(density < threhold);
index_be = min(index_set);
end