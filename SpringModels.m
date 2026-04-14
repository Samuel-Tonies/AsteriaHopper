
% Fuel slosh parameters
m_s = 355.9/2; % 50% fill
k_s = 500;
c_s = 2;

% Floor parameters
k_T = 9999999999;
c_T = 1500000;
mu_floor = .2;
% Define floor vector
floor_angle = 6; %deg
% n_floor = [cosd(90+floor_angle), 0, sind(90+floor_angle)];
n_floor = [0 0 1];
% Leg Model
angle_legs = 54; %degrees
body_length = 1.61544;
leg_length = 1.87198;

r_l1 = [cosd(angle_legs); 0; sind(angle_legs)];
r_l1 = r_l1/norm(r_l1);
r_l2 = [-cosd(angle_legs); 0; sind(angle_legs)];
r_l2 = r_l2/norm(r_l2);
r_l3 = [0; cosd(angle_legs); sind(angle_legs)];
r_l3 = r_l3/norm(r_l3);
r_l4 = [0; -cosd(angle_legs); sind(angle_legs)];
r_l4 = r_l4/norm(r_l4);

p_l1 = [-1; 0; -.5];
p_l2 = [1; 0; -.5];
p_l3 = [0; -1; -.5];
p_l4 = [0; 1; -.5]; %UPDATE ALL

k_l = 10000000000;
c_l = 2*sqrt(k_l*m_leg);
m_leg = 10;
leg_length = 1;