%% global data
mu = 398600;
%%%% part a %%%%
%% first orbit circular
r_a_1 = 6570; % circular
r_p_1 = 6570; % circular

h_1 = h_calculator(mu, r_a_1, r_p_1);

%% second orbit ellipse
r_p_2 = 6570 ;
r_a_2 = 42160;

h_2 = h_calculator(mu, r_a_2, r_p_2);


%% second orbit circular
r_a_3 = 42160; % circular
r_p_3 = 42160; % circular

h_3 = h_calculator(mu, r_a_3, r_p_3);

%% delta v

delta_v_1 = (h_2 - h_1) / r_a_2;

delta_v_2 = (h_3 - h_2) / r_p_2;

delta_v = delta_v_1 + delta_v_2;



%%%% part b %%%%
e = (r_p_2 - r_a_2) / (r_p_2 + r_a_2);

theta = pi;

a = (r_a_2 + r_p_2) /2;

tau = 2 * pi * sqrt(a^3 / mu);

time =  tau/2;








function h = h_calculator(mu, ra, rp)
    h = sqrt(2 * mu * (ra * rp) / (ra + rp));
end
