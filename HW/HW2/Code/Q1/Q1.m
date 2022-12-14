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
%% part c %%
v_mean = (v_1 + v_2_a) / 2;
distance = v_mean * delta_t;
delta_distance = delta_v * delta_t;
%% part d %%
tau_1 = sqrt(mean(r_1+6378)^3 / mu);
delta_t_1 = delta_t / tau_1;

tau_2 = sqrt(mean(r_2+6378)^3 / mu);
delta_t_2 = delta_t / tau_2;

%% part e %%
c_1 = 2*pi*(r_1(1)+6378);
delta_d = delta_distance / c_1;










%% functions
function h = h_calculator(mu, ra, rp)
    h = sqrt(2 * mu * (ra * rp) / (ra + rp));
end