%% global data
mu = 398600;
r = [1600, 5310, 3800]; %km
v = [-7.350, 0.4600, 2.470]; % km/s
%% part a %%
v_r = dot(r, v) / norm(r); % v_r > 0  away from perigee
h   = cross(r, v);
e   = (cross(v, h) - mu * r / norm(r))/mu;
theta = acos(dot(e, r) / norm(e) / norm(r));
a = norm(h)^2 / (mu * (1 - norm(e)^2));
tau = 2 * pi * sqrt(a^3 / mu);
M_e = E_calculator(theta, norm(e)) - ...
    norm(e) * sin(E_calculator(theta, norm(e)));

t = M_e / 2 / pi * tau;

%% tabulate

r = [1600, 5310, 3800]; %km
v = [-7.350, 0.4600, 2.470]; % km/s

[t_orbit_table, x_table] = ode45(@diff_eq_orbit,t:100:t+tau,[r';v']);

for i=1:length(t_orbit_table)
    fprintf("%.2f & %.2f & %.2f & %.2f  & %.2f & %.2f &  %.2f &  %.2f & %.2f \\\\ \n ",...
        t_orbit_table(i), x_table(i, 1), x_table(i, 2)...
        , x_table(i, 3), norm(x_table(i, 1:3)), x_table(i, 4),...
        x_table(i, 5), x_table(i, 6), norm(x_table(i, 4:6)))
end

%% better ode

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-9);
[t_orbit, x] = ode45(@diff_eq_orbit,[t, t+tau],[r';v'], options);

r = x(:, 1:3);
v = x(:, 4:end);

norm_r = zeros(length(x(:, 1)), 1);
for i = 1:length(norm_r)
    norm_r(i) = norm(r(i, :));
end

norm_v = zeros(length(x(:, 1)), 1);
for i = 1:length(norm_r)
    norm_v(i) = norm(v(i, :));
end


[xs, ys, zs] = sphere(50);

surf(6378*xs, 6378*ys, 6378*zs)
axis equal
colormap('gray');
xlabel('X(km)')
ylabel('Y(km)')
zlabel('Z(km)')

hold on
plot3(r(:, 1), r(:, 2), r(:, 3), 'b', 'linewidth', 4)

print('../../Figure/Short_project/3Dof_view.png','-dpng','-r400');

view(0,0)

print('../../Figure/Short_project/xz_view.png','-dpng','-r400');

view(-90,0)

print('../../Figure/Short_project/yz_view.png','-dpng','-r400');

view(-90,-90)

print('../../Figure/Short_project/xy_view.png','-dpng','-r400');



%% part b %%

N = cross([0, 0, 1], h);


if (e(3)) > 0
    omega = acos(dot(N / norm(N), e / norm(e)));
else
    omega = 2*pi - acos(dot(N / norm(N), e / norm(e)));
end



% [phi, lambda] = latlon(r);

for i=1:length(r)
    r(i, :) = r(i, :) * T_m;
end
phi = zeros(1, length(r));
lambda = zeros(1, length(r));
for i=1:length(r)
    [phi(i), lambda(i)] = latlon(r(i, :));
end














function E = E_calculator(theta, e)
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2));
end

function d = diff_eq_orbit(~, x)
mu = 398600;

r     = x(1:3);
r_dot = x(4:end);
d = zeros(6, 1);


d(1:3)   = r_dot;
d(4:end) = -mu / norm(r)^3 * r;
end


function [phi, lambda] = latlon(r)
    phi = asin(r(3)/norm(r));
    if r(2) > 0
        lambda = acos(r(1) / norm(r(1:2)));
    else
        lambda = 2*pi - acos(r(1) / norm(r(1:2)));
    end
end