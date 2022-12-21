%% ground track %%
%% part b %%
load("r_v_data.mat", 'r');
load("data_Q2_b.mat", "t");


for i=1:length(r)
    r(i, :) = r(i, :) * tarns_matrix(t(i));
end
phi = zeros(1, length(r));
lambda = zeros(1, length(r));
for i=1:length(r)
    [phi(i), lambda(i)] = latlon(r(i, :));
end

plot((lambda-pi)*180/pi, phi*180/pi, '.')
xlabel('longitude $\lambda^{\circ}$', 'Interpreter','latex',...
    'FontSize', 20);
ylabel('latitude $\phi^{\circ}$', 'Interpreter','latex', ...
    'FontSize', 20);
axis equal
axis([-180 180 -90 90])


print -depsc ../../Figure/Q2/latlong



%% on earth fig
% initializes figure
figure('Position',[300,300,1000,500]);
% plotting options
opts_earth1.Color = [140,21,21]/255;
opts_earth1.LineWidth = 2.5;
ground_track((phi)*180/pi,(lambda-pi)*180/pi,opts_earth1,'Earth');

print -depsc ../../Figure/Q2/latlong_earth


function [phi, lambda] = latlon(r)
    phi = asin(r(3)/norm(r));
    if r(2) > 0
        lambda = acos((r(1)/norm(r))/ cos(phi));
    else
        lambda = 2*pi - acos((r(1)/norm(r)) / cos(phi));
    end
end


function T_s =  tarns_matrix(t)
    omega = 2 * pi / 24 / 3600;
    T_s = [cos(omega*t), -sin(omega*t), 0;
           sin(omega*t),  cos(omega*t), 0;
                0              0      , 1];
end