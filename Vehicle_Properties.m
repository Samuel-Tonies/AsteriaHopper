m_l = 803; %kg
m_e = 150;
g = 1.625; 
l = 4.2926/2;

k_T = .8; % floor spring coefficient
C_T = 0.1; % floor damp coefficient
mu_surface = 0.2; % floor friction coefficient
n_surface = [0; 0; 1]; % surface plane normal vector

angle_legs = 54; %degrees
body_length = 1.61544;
leg_length = 1.87198;

r_l1 = [cosd(angle_legs); 0; sind(leg_length)];
r_l2 = [-cosd(angle_legs); 0; sind(leg_length)];
r_l3 = [0; cosd(angle_legs); sind(leg_length)];
r_l4 = [0; -cosd(angle_legs); sind(leg_length)];

