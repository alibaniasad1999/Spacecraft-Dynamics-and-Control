I_x = 8 ; %kg.m^2
I_y = 10; %kg.m^2
I_z = 14; %kg.m^2

sigma_x = (I_y - I_z) / I_x;
sigma_z = (I_y - I_x) / I_z;

omega_0 = sqrt(398500 / (6378+500)^3);

s = tf('s');

G = s^4+omega_0^2*(3*sigma_x+sigma_x*sigma_z+1)*s^2 ...
    +4*omega_0^2*sigma_x*sigma_z;

G = 1/G;

plot(pole(G), '*', 'MarkerSize',10, 'LineWidth',1)


set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
title('')
xlabel('Real Axis($\sec^{-1}$)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('Imaginary Axis($\sec^{-1}$)', 'interpreter', 'latex', 'FontSize', 24);
grid on
print('../../Figure/Q4/pole_zero_map','-depsc');



