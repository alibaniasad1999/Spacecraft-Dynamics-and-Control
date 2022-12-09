r_ap = [9798 5657 11314];
v_r = [-0.1 0.2 0.3];
mu = 398600;
a = norm(r_ap)/2;
theta = acos(dot(r_ap, v_r)/ norm(r_ap) / norm(v_r));

syms h e
eq1 = h^2 / mu == a*(1-e^2);
eq2 = norm(v_r) == mu / h * e * sin(theta);

answer_eq = solve([eq1, eq2]);

e = double(answer_eq.e(2));
h = double(answer_eq.h(2));

v_n = h / norm(r_ap);

gamma = atan(norm(v_r)/v_n);

r = h^2 / mu / (1+e*cos(theta));