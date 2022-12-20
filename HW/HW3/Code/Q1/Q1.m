clear;
clc;
%% Q1 %%

r = [0, 0, 0]; %km
v = [-0.1, -0.04, -0.02]; % km/s

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-9);
ode_time = 0:0.01:320*60; % sec
[t, x] = ode45(@diff_eq_orbit, ode_time, [r'; v'], options);

plot(t/60, x(:, 1:3), 'LineWidth',2)
xlim([0, 320]);
legend('x', 'y', 'z', 'Location','northwest', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Time (minute)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('Position (meter)', 'interpreter', 'latex', 'FontSize', 24);
print('../../Figure/Q1/position','-depsc');

plot(t/60, x(:, 4:end), 'LineWidth',2)
xlim([0, 320]);
ylim([-0.2 0.6])
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$', 'interpreter', 'latex', 'Location','northeast', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Time (minute)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity (meter/second)', 'interpreter', 'latex', 'FontSize', 24);
print('../../Figure/Q1/velocity','-depsc');






function d = diff_eq_orbit(~, x)

r_a = 590; % km
R_e = 6378; % km
r = r_a + R_e;
mu = 398600;
n = sqrt(mu / r^3);

d = zeros(6, 1);

d(1) = x(4); % dot x
d(2) = x(5); % dot y
d(3) = x(6); % dot z

d(4) = 2 * n * x(5) + 3 * n^2 * x(1); % ddot x
d(5) = -2 * n * x(4); % ddot y
d(6) = -n^2 * x(3); % ddot z

end