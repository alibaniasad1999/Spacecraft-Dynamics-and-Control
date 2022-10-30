clear;
clc;

%% global data
mu = 398600;

%% part a
ft2km = 0.0003048;

r = [4.1852, 6.2778, 10.463] * 10^7;
r = r * ft2km;

v = [2.5936, 0, 5.1872] * 10^4;
v = v * ft2km;

h = cross(r, v);

e = cross(v, h) / mu - r / norm(r);

epsilon = norm(v)^2 / 2 - mu / norm(r);

gamma = asin(dot(r, v) / norm(r) / norm(v));


fprintf("gamma(deg) = %f \n", gamma * 180 / pi);

P = norm(h)^2 / mu;


