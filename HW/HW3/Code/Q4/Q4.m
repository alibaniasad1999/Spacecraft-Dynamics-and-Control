%% part a %%
G = 6.67259e-20;
m_1 = 5.974e24;
m_2 = 7.348e22;
mu_1 = G * m_1;
mu_2 = G * m_2;

r_12 = 3.844e5;

omega = sqrt(G * (m_1 + m_2) / r_12^3);

r = [0.994, 0, 0];
v = [0, -2.001585106, 0];
mu_u = m_2 / (m_1 + m_2);
C = 1/2 * norm(v)^2 - 1/2*(r(1)^2 + r(2)^2) - mu_u;






%% canonical %%

global mu r_1 r_2



mu = m_2 / (m_1 + m_2);

r_1 =  -[mu, 0, 0] + r;
r_2 = [1-mu, 0, 0] + r;

r_1 = norm(r_1);
r_2 = norm(r_2);


[t, X] = ode45(@diff_eq_orbit_global_canonical,0:0.001:100,[r';v']);

%% ploter %%


plot(mu, 0, '.', 'markersize', 160)
hold on
plot(mu-1, 0, '.', 'markersize', 2)
plot(X(:, 1), X(:, 2))

text(-0.05, 0, 'Earth')
xlabel('X', 'interpreter', 'latex', 'FontSize', 24);
ylabel('Y', 'interpreter', 'latex', 'FontSize', 24);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
print('../../Figure/Q4/xyplane','-depsc');


%% part c %%
X_0 = 0.994;
Y_0 = 0;
U_xx = 1 - ((1-mu)*(1/r_1^3-3*(X_0 - mu)^2/r_1^5)+mu*(1/r_2^3-3*(X_0+1-mu)^2/r_2^5));
U_yy = 1 - ((1-mu)*(1/r_1^3-3*Y_0^2/r_1^5)+ mu*(1/r_2^3-3*Y_0^2/r_2^5));

syms lambda
eq_lambda = lambda^4 + (4-U_xx-U_yy)\lambda^2 + U_xx*U_yy == 0; 
% lambda_ans = solve(eq_lambda);
lambda_ans = double(vpa(solve(eq_lambda)))
%% subfunctions %%
function d = diff_eq_orbit_global_canonical(~, x)

global mu r_1 r_2

d = zeros(6, 1);

d(1) = x(4); % dot x
d(2) = x(5); % dot y
d(3) = x(6); % dot z

d(4) = 2*x(5) + x(1) - (1 - mu)/r_1^3 * (x(1)-mu) - mu/r_2^3*(x(1)+1-mu); % ddot x
d(5) = -2*x(4) + x(2) - (1-mu)/r_1^3* x(2) - mu/r_2^3 * x(2); % ddot y
d(6) = -(1-mu)/r_1^3 * x(3) - mu/r_2^3*x(3);
end