clc;
clear;
%% data %%
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

[V, ~] = eig(I); % eigenvectors 
[yaw, pitch, roll] = dcm2angle(V);

fprintf('roll is: %.2f\n', roll*180/pi);
fprintf('pitch is: %.2f\n', pitch*180/pi)
fprintf('yaw is: %.2f\n', yaw*180/pi)

I_principal_axes = V' * I * V;

omega = V' * omega;

h = I_principal_axes * omega;

T =  1/2 * omega' * I_principal_axes * omega;

%% Bonus %%
a = double(vpa(norm(h) / I_principal_axes(1, 1)));
b = double(vpa(norm(h) / I_principal_axes(2, 2)));
c = double(vpa(norm(h) / I_principal_axes(3, 3)));
alpha = sqrt(double(vpa(2*T / I_principal_axes(1, 1))));
beta = sqrt(double(vpa(2*T / I_principal_axes(2, 2))));
gamma = sqrt(double(vpa(2*T / I_principal_axes(3, 3))));

