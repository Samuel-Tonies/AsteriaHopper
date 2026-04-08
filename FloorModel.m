% Trying to model floor

% Define floor vector
floor_angle = 6; %deg
n_floor = [cosd(90+floor_angle), 0, sind(90+floor_angle)];


angle_legs = 54; %degrees
body_length = 1.61544;
leg_length = 1.87198;

r_l1 = [cosd(angle_legs); 0; sind(leg_length)];
r_l2 = [-cosd(angle_legs); 0; sind(leg_length)];
r_l3 = [0; cosd(angle_legs); sind(leg_length)];
r_l4 = [0; -cosd(angle_legs); sind(leg_length)];