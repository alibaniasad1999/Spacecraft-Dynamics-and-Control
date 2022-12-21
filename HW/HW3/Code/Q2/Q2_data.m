%% Q2 data %%
r = [1600, 5310, 3800]; % km
v = [-7.35, 0.46, 2.47]; % km/s
h = cross(r, v);
v_r = dot(r, v) / norm(r);
mu = 398600;
e = (cross(v, h) - mu * r / norm(r)) / mu;

a = norm(h)^2 / (mu * (1 - norm(e)^2));

N = cross([0, 0, 1], h);


%% theta %%
if v_r >= 0
    theta = acos(dot(e, r) / norm(e) / norm(r));
else
    theta = 2*pi - acos(dot(e, r) / norm(e) / norm(r));
end
%% Omega %%
if N(2) >= 0
    Omega = acos(N(1)/norm(N));
else
    Omega = 2*pi - acos(N(1)/norm(N));
end

%% omega %%
if e(3) >= 0
    omega = acos(dot(N, e)/norm(e)/norm(N));
else
    omega = 2*pi - acos(dot(N, e)/norm(e)/norm(N));
end
%% i %%
i = acos(h(3)/norm(h));
%% use function %%
Q2(a, norm(e), i, Omega, theta)
% (a_data, e_data, i_data, Omega_data, theta_data)