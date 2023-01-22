clc;
clear;
%% Q1 part I %%
syms I_yz
I = [30,  -10,    0;
    -10,   20,-I_yz;
      0,-I_yz,  30];
omega = 10 * ones(3, 1);
h = [200, 200, 400]';

ans_I_yz = solve(h == I * omega);

I = subs(I, I_yz, double(ans_I_yz));

I = double(I);

T =  1/2 * omega' * I * omega;

fprintf('RKE is: %.2f\n', T)

I_xi = 2 * T / dot(omega,omega);

fprintf('I_xi is: %.2f\n', I_xi)
%% part c %%

[V, ~] = eig(I); % eigenvectors 
[yaw, pitch, roll] = dcm2angle(V);

fprintf('roll is: %.2f\n', roll*180/pi);
fprintf('pitch is: %.2f\n', pitch*180/pi)
fprintf('yaw is: %.2f\n', yaw*180/pi)

I_principal_axes = V' * I * V;

omega = V' * omega;

h = I_principal_axes * omega;

T =  1/2 * omega' * I_principal_axes * omega;

%% elipsoid if inertia %%

[x, y, z] = ...
ellipsoid(0, 0, 0, sqrt(1 / I_principal_axes(1, 1)),...
                   sqrt(1 / I_principal_axes(2, 2)),...
                   sqrt(1 / I_principal_axes(3, 3)),50);

colormap('winter');
surf(x, y, z);
shading interp
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('X', 'interpreter', 'latex', 'FontSize', 24);
ylabel('Y', 'interpreter', 'latex', 'FontSize', 24);
zlabel('Z', 'interpreter', 'latex', 'FontSize', 24);


print -depsc ../../Figure/Q1/3Dof_view_elipsoid_inertia

view(0,0)

print -depsc ../../Figure/Q1/xz_view_elipsoid_inertia

view(-90,0)

print -depsc ../../Figure/Q1/yz_view_elipsoid_inertia

view(-90,-90)

print -depsc ../../Figure/Q1/xy_view_elipsoid_inertia



%% ploter %% 
w = winter;

[x, y, z] = ...
ellipsoid(0, 0, 0, double(vpa(norm(h) / I_principal_axes(1, 1))),...
                   double(vpa(norm(h) / I_principal_axes(2, 2))),...
                   double(vpa(norm(h) / I_principal_axes(3, 3))),50);

colormap('autumn');
surf(x, y, z);

hold on

[x, y, z] = ...
ellipsoid(0, 0, 0, sqrt(double(vpa(2*T / I_principal_axes(1, 1)))),...
                   sqrt(double(vpa(2*T / I_principal_axes(2, 2)))),...
                   sqrt(double(vpa(2*T / I_principal_axes(3, 3)))),1000);

plot3(x, y, z, 'color' , "#0072BD", 'LineWidth', 10);
% colormap('autumn');
axis equal

shading interp

print -depsc ../../Figure/Q1/3Dof_view

view(0,0)

print -depsc ../../Figure/Q1/xz_view

view(-90,0)

print -depsc ../../Figure/Q1/yz_view

view(-90,-90)

print -depsc ../../Figure/Q1/xy_view

close all




