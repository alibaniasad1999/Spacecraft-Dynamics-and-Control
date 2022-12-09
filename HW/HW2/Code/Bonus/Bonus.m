syms theta1 theta2
eq1 = 7200 / (1+0.5*cos(theta1)) == 8064 / (1+0.4*cos(theta2));
theta1_sol = solve(eq1, theta1);

v1 = [-sin(theta1_sol), cos(theta1_sol)];
v2 = [-sin(theta2), cos(theta2)];

func = (v1 - v2).^2;

cost_1 = @(theta2)(sin(theta2) - (1 - ((5*cos(theta2))/7 - 3/14)^2)^(1/2))^2 +...
    ((2*cos(theta2))/7 + 3/14)^2;


cost_2 = @(theta2) + (sin(theta2) + (1 - ((5*cos(theta2))/7 - 3/14)^2)^(1/2))^2 +...
    ((2*cos(theta2))/7 + 3/14)^2;
x = fminbnd(cost_2,-2*pi,2*pi);
data = ones(1, length(dim));
for i = 1:length(dim)
data(i) = cost_1(dim(i));
end
plot(dim , data)
