%% ploter part b%%
load('data_Q2.mat')
r = zeros(length(y), 3);
v = zeros(length(y), 3);
mu = 398600;
for i = 1:length(y)
    [r(i, :), v(i, :)] = sv_from_coe(y(i, :), mu);
end
save('r_v_data_all', "v", "r");

load('data_Q2_b.mat')
r = zeros(length(y), 3);
v = zeros(length(y), 3);
mu = 398600;
for i = 1:length(y)
    [r(i, :), v(i, :)] = sv_from_coe(y(i, :), mu);
end
save('r_v_data', "v", "r");