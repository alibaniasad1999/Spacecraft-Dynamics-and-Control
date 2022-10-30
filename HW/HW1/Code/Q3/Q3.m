clear;
clc;

%% global data
R_e = 6378; 
mu = 398600;
%% Q3 part a %%
r_p = R_e + 500;
r_a = R_e + 5000;
e = (r_a - r_p) / (r_a + r_p);
a = (r_a + r_p) / 2;

tau = 2 * pi * sqrt(a^3 / mu);

%% solve theta
theta = sym('theta','real');
eq = e * cos(theta) - (1 - e^2) * a / R_e * sin(theta) + 1;
theta_ans = solve(eq, theta);
theta_ans_vpa = vpasolve(eq, theta);
theta = [vpa(theta_ans(1)), theta_ans_vpa];
theta = double(theta);
M_e_1 = E_calculator(theta(1), e) - e * sin(E_calculator(theta(1), e));
M_e_2 = E_calculator(theta(2), e) - e * sin(E_calculator(theta(2), e));

t_1 = M_e_1 / 2 / pi * tau;
t_2 = M_e_2 / 2 / pi * tau;

fprintf("time part a: %f \n", abs(t_1 - t_2));

%% Q3 part a %%
theta = sym('theta','real');
eq = e * cos(theta) - (1 - e^2) * a / R_e * cos(theta) + 1;
theta_ans_1 = solve(eq, theta);
eq = e * cos(theta) + (1 - e^2) * a / R_e * cos(theta) + 1;
theta_ans_2 = solve(eq, theta);
theta = double(vpa([theta_ans_1(2), theta_ans_2(2)]));

M_e_12 = E_calculator(theta(1), e) - e * sin(E_calculator(theta(1), e));
M_e_22 = E_calculator(theta(2), e) - e * sin(E_calculator(theta(2), e));

t_12 = M_e_12 / 2 / pi * tau;
t_22 = M_e_22 / 2 / pi * tau;

fprintf("time part b: %f \n", abs(t_12 - t_22));






function E = E_calculator(theta, e)
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2));
end