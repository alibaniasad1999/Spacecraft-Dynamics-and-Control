
% global omega mu_1 mu_2 r_1 r_2 r_12 pi_1 pi_2


% G = 6.6742e-11;
% 
% m_1 = 5.974e24;
% m_2 = 7.348e22;
% r_12 = 3.844e5;
% 
% r_1 = norm([0.994, 0, 0])*r_12;
% r_2 = norm([0.994, -2.001585106, 0])*r_12;
% 
% pi_1 = m_1 / (m_1 + m_2);
% pi_2 = m_2 / (m_1 + m_2);
% 
% mu = G*(m_1 + m_2);
% mu_1 = G * m_1;
% mu_2 = G * m_2;
% 
% omega = sqrt(mu / r_12^3);



% [t, X] = ode45(@diff_eq_orbit,0:10000,[r';v']);








% function d = diff_eq_orbit_global(~, x)
% global omega mu_1 mu_2 r_1 r_2 r_12 pi_1 pi_2
% 
% d = zeros(6, 1);
% 
% d(1) = x(4); % dot x
% d(2) = x(5); % dot y
% d(3) = x(6); % dot z
% 
% d(4) = 2*omega*x(5) + omega^2*x(1) - mu_1/r_1^3*(x(1)-pi_2*r_12) - ...
%       -mu_2/r_2^3*(x(1)-pi_1*r_12); % ddot x
% d(5) = -2*omega*x(4) + omega^2*x(2) - mu_1/r_1^3*x(2) - mu_2/r_2^3*x(2); % ddot y
% d(6) = -mu_1/r_1^3*x(3) - mu_2/r_2^3*x(3);
% end




%% canonical %%

global mu r_1 r_2

m_1 = 5.974e24;
m_2 = 7.348e22;
r_12 = 3.844e5;

r = [0.994, 0, 0];
v = [0, -2.001585106, 0];

mu = m_2 / (m_1 + m_2);

r_1 =  -[mu, 0, 0] + r;
r_2 = [1-mu, 0, 0] + r;

r_1 = norm(r_1);
r_2 = norm(r_2);


[t, X] = ode45(@diff_eq_orbit_global_canonical,0:0.001:100,[r';v']);


function d = diff_eq_orbit_global_canonical(~, x)

global mu r_1 r_2

d = zeros(6, 1);

d(1) = x(4); % dot x
d(2) = x(5); % dot y
d(3) = x(6); % dot z

d(4) = 2*x(5) - x(1) - (1 - mu)/r_1^3 * (x(1)-mu) - mu/r_2^3*(x(1)+1-mu); % ddot x
d(5) = -2*x(4) + x(2) - (1-mu)/r_1^3* x(2) - mu/r_2^3 * x(2); % ddot y
d(6) = -(1-mu)/r_1^3 * x(3) - mu/r_2^3*x(3);
end