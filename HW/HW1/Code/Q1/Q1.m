%%%%%%%%%% Q1 %%%%%%%%%% 

%% a
r = [0 2 0];
v = 1 / sqrt(2) * [ 1 1 0];

h = cross(r, v);

C = cross(v, h) - r / norm(r);

mu = 1;

e = C / mu;

fprintf("h.e = %d \n", floor(dot(h, e)))

theta = acos(norm(h)^2 / mu / norm(r) - 1) / norm(e);

fprintf("Theta(deg) = %f \n", theta * 180 / pi);


epsilon = norm(v)^2 / 2 - mu / 2;

r_2 = 32;
v_2 = sqrt(2 * mu / r_2);

theta_2 = acos(norm(h)^2 / mu / r_2 - 1) / norm(e);

fprintf("Theta2(deg) = %f \n", theta_2 * 180 / pi);

