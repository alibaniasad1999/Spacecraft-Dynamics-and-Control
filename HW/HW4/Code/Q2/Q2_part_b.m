clc;
clear;
%% part a %%
q = [0.3, 0.2, 0.5, 0.7874];
euler_angles = quat2eul(q) * 180/pi;

roll = euler_angles(1);
pitch = euler_angles(2);
yaw = euler_angles(3);

fprintf('roll: %.2f\n', roll);
fprintf('pitch: %.2f\n', pitch);
fprintf('yaw: %.2f\n', yaw);

%% part c %%

omega = [2,3,5] * 10^-3 * 2 * pi;

initial_condition = zeros(6, 1);

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-9);
[t_non_linear, x_non_linear] = ode45(@diff_eq_non_linear,0:50,initial_condition, options);


[t_linear, x_linear] = ode45(@diff_eq_linear,0:50,initial_condition, options);

%% ploter %%
% phi part_d
plot(t_linear, x_linear(:, 4)*180/pi, 'LineWidth',2)
hold on
plot(t_non_linear, x_non_linear(:, 4)*180/pi, 'LineWidth',2)
legend('linear', 'non linear', 'Location','northwest', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Time ($\sec$)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\phi^{\circ}$', 'interpreter', 'latex', 'FontSize', 24);
print('../../Figure/Q2/phi_part_d','-depsc');
hold off
% theta part_d
plot(t_linear, x_linear(:, 5)*180/pi, 'LineWidth',2)
hold on
plot(t_non_linear, x_non_linear(:, 5)*180/pi, 'LineWidth',2)
legend('linear', 'non linear', 'Location','northwest', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Time ($\sec$)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\theta^{\circ}$', 'interpreter', 'latex', 'FontSize', 24);
print('../../Figure/Q2/theta_part_d','-depsc');
hold off
% psi part_d
plot(t_linear, x_linear(:, 6)*180/pi, 'LineWidth',2)
hold on
plot(t_non_linear, x_non_linear(:, 6)*180/pi, 'LineWidth',2)
legend('linear', 'non linear', 'Location','northwest', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Time ($\sec$)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\psi^{\circ}$', 'interpreter', 'latex', 'FontSize', 24);
print('../../Figure/Q2/psi_part_d','-depsc');
hold off
%% ode functions %%
function d = diff_eq_non_linear(~, x)

% data needed %
I_x = 8 ;
I_y = 10;
I_z = 14;

R_0 = 6378 + 500;
mu  = 398500;
omega_0 = sqrt(mu/R_0^3);


omega_x = x(1);
omega_y = x(2);
omega_z = x(3);

phi   = x(4);
theta = x(5);
psi   = x(6);

G = 3*mu/2/R_0^3*[(I_z-I_y)*sin(2*phi)*cos(theta)^2;...
                  (I_z-I_x)*sin(2*theta)*cos(phi) ;...
                  (I_x-I_y)*sin(2*theta)*sin(phi)];

G_x = G(1);
G_y = G(2);
G_z = G(3);

M_x = 8 * 10^-6;
M_y = 8 * 10^-6;
M_z = 8 * 10^-6;

d = zeros(6, 1);



d(1) = (G_x + M_x + omega_y * omega_z * (I_y - I_z)) / I_x;
d(2) = (G_y + M_y + omega_x * omega_z * (I_z - I_x)) / I_y;
d(3) = (G_z + M_z + omega_y * omega_x * (I_x - I_y)) / I_z;

C_r2b = [cos(theta)*cos(psi) , -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi) , sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi); ... 
        cos(theta)*sin(psi) , cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi) , -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi); ...
        -sin(theta) , sin(phi)*cos(theta) , cos(phi)*cos(theta)]';

pqr = d(1:3) + C_r2b * [0; omega_0; 0];

p = pqr(1);
q = pqr(2);
r = pqr(3);

euler_propagation = [1, sin(phi)*tan(theta), cos(theta)*tan(theta);...
                    0,      cos(theta),          -sin(theta);...
                    0,  sin(phi)/cos(theta), cos(phi)/cos(theta)];

d(4:6) = euler_propagation * [p; q; r];

end


function d = diff_eq_linear(~, x)

% data needed %
I_x = 8 ;
I_y = 10;
I_z = 14;

R_0 = 6378 + 500;
mu  = 398500;
omega_0 = sqrt(mu/R_0^3);


omega_x = x(1);
omega_y = x(2);
omega_z = x(3);

phi   = x(4);
theta = x(5);
psi   = x(6);

G = 3*mu/2/R_0^3*[(I_z-I_y)*2*phi;...
                  (I_z-I_x)*2*theta;...
                  0];

G_x = G(1);
G_y = G(2);
G_z = G(3);

M_x = 8 * 10^-6;
M_y = 8 * 10^-6;
M_z = 8 * 10^-6;

d = zeros(6, 1);



d(1) = (G_x + M_x + omega_y * omega_z * (I_y - I_z)) / I_x;
d(2) = (G_y + M_y + omega_x * omega_z * (I_z - I_x)) / I_y;
d(3) = (G_z + M_z + omega_y * omega_x * (I_x - I_y)) / I_z;

C_r2b = [1, -psi , theta; ... 
        psi , 1, -phi; ...
        -theta, phi,1]';

pqr = d(1:3) + C_r2b * [0; omega_0; 0];

d(4:6) = pqr;

end
