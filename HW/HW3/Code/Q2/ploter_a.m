%% ploter a %%
function ploter_a(a_data, e_data, i_data, Omega_data, theta_data)
load('data_Q2.mat', 'y', 't');
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians
%...Constants;
mu = 398600; %Gravitational parameter (km^3/s^2)
% RE = 6378; %Eath's radius (km)
% c = 2.998e8; %Speed of light (m/s)
% S = 1367; %Solar constant (W/m^2)
% Psr = S/c; %Solar pressure (Pa);

%...Satellite data:
% CR = 2; %Radiation pressure codfficient
% m = 100; %Mass (kg)
% As = 200; %Frontal area (m^2);
%...Initial orbital parameters (given):
a0 = a_data; %Semimajor axis (km)
e0 = e_data; %eccentricity
incl0 = i_data; %Inclination (radians)
RA0 = Omega_data; %Right ascencion of the node (radians)
TA0 = theta_data; %True anomaly (radians)
w0 = 227.493*deg; %Argument of perigee (radians)
%...Initial orbital parameters (inferred):
h0 = sqrt(mu*a0*(1-e0^2)); %angular momentrum (km^2/s)
% T0 = 2*pi/sqrt(mu)*a0^1.5; %Period (s)
% rp0 = h0^2/mu/(1 + e0); %perigee radius (km)
% ra0 = h0^2/mu/(1 - e0); %apogee radius (km)
%...Store initial orbital elements (from above) in the vector coe0:

h = y(:,1);
e = y(:,2);
RA = y(:,3);
incl = y(:,4);
w = y(:,5);
TA = y(:,6);

plot(t/days,h - h0, 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$h-h_0$', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/h_fig','-depsc');


plot(t/days,e - e0, 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$e-e_0$', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/e_fig','-depsc');


plot(t/days,(RA - RA0)/deg, 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\Omega-Omega_0$ (deg)', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/Omega_fig','-depsc');


plot(t/days,(incl - incl0)/deg, 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$i-i_0$ (deg)', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/i_fig','-depsc');


plot(t/days,(w - w0)/deg, 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\omega-\omega_0$ (deg)', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/omega_fig','-depsc');


plot(t/days,mod((TA - TA0)/deg, 360), 'LineWidth',2)
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Days', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\theta-\theta_0$ (deg)', 'interpreter', 'latex', 'FontSize', 24);
axis tight
print('../../Figure/Q2/theta_fig','-depsc');
end