%% Trajectory Analysis Tool
clc;
% Tool for analyzing and sizing hopper system
% Note: program may occasionally fail because of numerical solver issues,
% reach out and I can help

% Currently, there is no functionality to use real terrain for generating
% trajectories. As a result, there is a chance of accidentaly creating a
% trajectory that would clip walls or rim. So use with that knowledge

% Enter Desired Distances of Trajectory
x_des_m = 10000; % Horizontal distance from initial point in meters
z_des_m = -4000; % Vertical distance from initial point in meters

% Enter Hopper Properties
m0_kg = 900; % Total mass of hopper in kg
Isp_s = 220; % Isp of propellant system

% System Bounds
% It may be desireable to set limits on the maximum thrust or vehicle pitch
% Adjust system bounds as desired

% It is recommended to keep thrust and pitch bounded and just set
% reasonable values, but play with it as you desire

thrust_bound_b = true; % Disables or enables thrust limits
thrust_max_N = 9999; % Maximum thrust in Newtons
thrust_min_N = 0; % Minimum thrust in Newtons

% Theta is effectively the vehicle pitch
% I've read 9 degrees is a reasonable assumption for vertical takeoff
% vertical lander hopper
theta_bound_b = true;
theta_max_deg = 9;  % Maximum pitch in degrees

% Chose to represent ground or not in trajectory generation
% Right now ground is represented by sinusoidal wave that tries to mimick 
% shackelton, go to ground_alt() to play around with it

% May break results when enables
ground_on_b = false;

% Display Option
% Gives you option to see solver iteration if you so desire
show_solver_b = false;

%=========================== DO NOT EDIT ================================%
% Splice count
N = 40;

% Pack initial conditions
BC.x0 = 0;
BC.xdot0 = 0;
BC.z0 = 0;
BC.zdot0 = 0;

% Pack final conditions
BC.xf = x_des_m;
BC.xdotf = 0;
BC.zf = z_des_m;
BC.zdotf = 0;

% Pack Parameters
param.m = m0_kg;
param.g = 1.625;
if thrust_bound_b
    param.Fmax = thrust_max_N;
    param.Fmin = thrust_min_N;
else
    param.Fmax = 99*param.m;
    param.Fmin = 0;
end
if theta_bound_b
    param.theta_max = deg2rad(theta_max_deg);
    param.theta_min = -param.theta_max;
else
    param.theta_max = pi/2;
    param.theta_min = pi/2;
end
param.isp = 220;
param.ground = ground_on_b;
param.show_solver = show_solver_b;

fprintf('Thinking....\n')
[sol, exitflag] = Trap2DTraj(N, BC, param);
if exitflag == 1
    fprintf('Solved!\n')
    % Populate workspace
    x = sol(1:N+1);
    xdot = sol(N+2:2*(N+1));
    z = sol(2*(N+1)+1:3*(N+1));
    zdot = sol(3*(N+1)+1:4*(N+1));
    F = sol(4*(N+1)+1:5*(N+1));
    theta = sol(5*(N+1)+1:6*(N+1));
    mass_kg = sol(6*(N+1)+1:7*(N+1));
    tf = sol(end);
    
    prop_burned = param.m - mass_kg(end);
    fprintf('Total time (s) of trajectory = %.2f\n',tf)
    fprintf('Propellant mass (kg) used in hop = %.2f\n', prop_burned)
    
    trajPlot(sol,N,param)
else
    fprintf("Solver could not find optimal solution :(\n")
end