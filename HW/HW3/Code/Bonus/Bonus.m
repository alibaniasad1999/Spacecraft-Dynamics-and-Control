r_a = 590; % km
R_e = 6378; % km
r = r_a + R_e;
mu = 398600;
n = sqrt(mu / r^3);
t = 8 * 3600; % sec
r_0 = [143, 137.279, 17.766]';
v_0 = [-0.1, 0.27, -0.01]';
Phi_rr = [4-3*cos(n*t)  0  0;...
          6*(sin(n*t)-n*t)  1  0; ...
           0  0  cos(n*t)];



Phi_rv = [1/n*sin(n*t), 2/n*(1-cos(n*t)),  0;...
    2/n*(cos(n*t)-1) , 1/n*(4*sin(n*t)-3*n*t), 0;...
    0 , 0 , 1\n*sin(n*t)];

Phi_vr = [3*n*sin(n*t) , 2*n*sin(n*t) , 0;...
    6*n*(cos(n*t) -1) , 0 , 0;...
    0 , 0 , -n\sin(n*t)];

Phi_vv = [cos(n*t) , 2*sin(n*t) , 0;...
    -2*n*sin(n*t) , 4*n*cos(n*t)-3 , 0;...
    0 , 0 , cos(n*t)];



delta_v_plus = -Phi_rv^-1 * Phi_rr * r_0;
delta_v_f_minus = Phi_vr * r_0 + Phi_vv * (-Phi_rv^-1 * Phi_rr * r_0);

%% LQR method %%
global A B Q R n_d R_inv
A = [Phi_rr, zeros(3, 3); zeros(3, 3), Phi_vr];
B = [Phi_rv, zeros(3, 3); zeros(3, 3), Phi_vv];


Q	= 1*eye(6);
R	= 10*eye(6);
R_inv = R^-1;
H	= 40*eye(6);
tf	= 4*3600;

K0	= H;
n_d	= 6;
k0	= reshape(K0,n_d^2,1);

global K
[~,LQR_gain,~] = icare(A,B,Q,R);
K = LQR_gain;
% [t_K,K_arr] = ode45(@diff_eq_Riccati,linspace(tf, 1, 10000),k0);


x0	= [1 1 1 1 1 1]';
[t,x] = ode45(@diff_eq_states,[0,tf],x0);

figure(101)
hold on
plot(t,x)





%%

function d = diff_eq_states(~,x)
global A B R_inv K
u	= -R_inv*B'*K*x;
d	= A*x + B*u;
end



        
