clc;
clear;
%% initial data %%
mass = 12.5e3; % kg
mu = 398600;
R_e = 6378; % km
trust_simultaneously = 53400;
%% part a %%
r_1 = [300, 300];
r_a_1 = R_e + r_1(1);
r_p_1 = R_e + r_1(2);

h_1 = h_calculator(mu, r_a_1, r_p_1);
v_1 = h_1 / r_a_1;

r_2 = [250, 300];
r_p_2 = R_e + r_2(1);
r_a_2 = R_e + r_2(2);

h_2 = h_calculator(mu, r_a_2, r_p_2);
v_2_a = h_2 / r_a_2;

delta_v = abs(v_2_a - v_1);

%% part b %%
delta_t = delta_v / trust_simultaneously * mass;










%% functions
function h = h_calculator(mu, ra, rp)
    h = sqrt(2 * mu * (ra * rp) / (ra + rp));
end