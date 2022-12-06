clear;
clc;

%% global data
R_e = 6378; 
mu = 398600;
%% Q3 part a %%
r_p = 8000;
r_a = 13000;
e = (r_a - r_p) / (r_a + r_p);
a = (r_a + r_p) / 2;

tau = 2 * pi * sqrt(a^3 / mu);

%% solve theta
theta_1 = pi/6;
theta_2 = pi/2;

M_e_1 = E_calculator(theta_1, e) - e * sin(E_calculator(theta_1, e));
M_e_2 = E_calculator(theta_2, e) - e * sin(E_calculator(theta_2, e));

t_1 = M_e_1 / 2 / pi * tau;
t_2 = M_e_2 / 2 / pi * tau;

fprintf("time part a: %f sec \n", abs(t_1 - t_2));
delta_t = abs(t_1 - t_2);

h = h_calculator(mu, r_a, r_p);

R_1 = h^2 / mu * (1 / (1 + e * cos(theta_1))); % norm R
R_2 = h^2 / mu * (1 / (1 + e * cos(theta_2))); % norm R

r_1 = [R_1 * cos(theta_1), R_1 * sin(theta_1), 0];
r_2 = [R_2 * cos(theta_2), R_2 * sin(theta_2), 0];

string = 'pro';
[v1, v2] = lambert(r_1, r_2, delta_t, string);

v1_o = mu/h * [-sin(theta_1), (e+cos(theta_1)), 0];
v2_o = mu/h * [-sin(theta_2), (e+cos(theta_2)), 0];



function E = E_calculator(theta, e)
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2));
end

function h = h_calculator(mu, ra, rp)
    h = sqrt(2 * mu * (ra * rp) / (ra + rp));
end