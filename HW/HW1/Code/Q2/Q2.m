clear;
clc;

R_e = 6378;
h = 600;
mu = 398600;
r = R_e + h;
v_r = 3.5; % dot(r)
v_h = 7;
v = sqrt(3.5^2 + 7^2);
h = r * v_h;

epsilon = v^2 / 2 - mu / r;

e = sqrt(1 + 2 * epsilon * h^2/(mu^2));

P = h^2 / mu;

theta = acos((P-r)/r/e);





