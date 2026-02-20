% Trajectory optimizer using trapezoidal integration
clear;clc;
T_max = 99999;
T_min = 100;
theta_max = 0.15;
theta_min = -theta_max;

N = 40;

param.m = 1000;
param.g = 1.62;
param.d = 1;
param.I = 6000;
param.Fmax = T_max;
param.Fmin = 100;
param.theta_max = theta_max;
param.theta_min = theta_min;
param.isp = 220;

BC.x0 = 0;
BC.xdot0 = 0;
BC.z0 = groundAlt(BC.x0);
BC.zdot0 = 0;
% BC.phi0 = 0;
% BC.phidot0 = 0;

BC.xf = 2700;
BC.xdotf = 0;
BC.zf = groundAlt(BC.xf);
BC.zdotf = 0;
% BC.phif = 0;
% BC.phidotf = 0;

DV0 = init_guess(N,param,BC);
[c0, ceq0] = nonlcon_trap(DV0, N, param, BC);
fprintf('max inequality violation = %.3e\n', max(c0));
fprintf('max equality violation = %.3e\n', max(abs(ceq0)));


options = optimoptions('fmincon',...
  'Algorithm','sqp', ...
  'MaxFunctionEvaluations',1e6, ...
  'MaxIterations',2e4, ...
  'Display','iter', ...
  'FiniteDifferenceType','central');  % central diffs are more accurate (costly)


sol = fmincon(@(DV) obj(DV,N),DV0, [], [], [], [], [], [], ...
                       @(DV)nonlcon_trap(DV,N,param,BC), options);

x = sol(1:N+1);
xd = sol(N+2:2*(N+1));
z = sol(2*(N+1)+1:3*(N+1));
zd = sol(3*(N+1)+1:4*(N+1));
F = sol(4*(N+1)+1:5*(N+1));
theta = sol(5*(N+1)+1:6*(N+1));
mass = sol(6*(N+1)+1:7*(N+1));
h = groundAlt(x);

prop_burned = param.m - mass(end);
fprintf('Propellant mass (kg) used in hop = %.2f\n', prop_burned)

function DV0 = init_guess(N,param,BC)
g = param.g;
xf = BC.xf;
tf0 = sqrt(2*xf/g); %1.5 times ballistic TOF
V_init = sqrt(xf*g);
U = param.isp*9.81;
prop_guess = param.m*(1-exp(-2*V_init/U));

t_guess = linspace(0,tf0,N+1)';
x_init = V_init*cosd(45)*t_guess;
xd_init = ones(N+1,1)*V_init*cosd(45);
z_init = g/2.*t_guess.^2+V_init*sind(45).*t_guess;
zd_init = g.*t_guess + V_init*sind(45);

% phi_init = zeros(N+1,1);
% phid_init = zeros(N+1,1);

F_init = zeros(N+1,1);
for i=1:ceil((N+1)/2)
        F_init(i) = param.Fmax*(1-2*t_guess(i)/tf0);
end
for i=ceil((N+1)/2)+1:N+1
        F_init(i) = F_init(i-ceil((N+1)/2));
end
    
theta_init = atan2(zd_init,xd_init);
m_init = param.m - t_guess/tf0*prop_guess;

% DV0 = [x_init;xd_init; z_init; zd_init; phi_init; phid_init; F_init;...
%     theta_init; tf0];
DV0 = [x_init;xd_init; z_init; zd_init; F_init;...
    theta_init; m_init; tf0];
end

function J = obj(DV,N)
    nNodes = N+1;
    F = DV(4*nNodes+1:5*nNodes);
    tf = DV(end);
    dt = tf / N;
    J = 0.5 * dt * sum(F(1:end-1).^2 + F(2:end).^2);
end

function [c, ceq] = nonlcon_trap(DV,N,param,BC)
% NONLCON_TRAP  nonlinear constraints (trapezoidal collocation) for 2D hopper
% --- Unpack states/controls/time ---
nNodes = N + 1;
x       = DV(1:nNodes);
x_dot   = DV(nNodes+1:2*nNodes);
z       = DV(2*nNodes+1:3*nNodes);
z_dot   = DV(3*nNodes+1:4*nNodes);
% phi     = DV(4*nNodes+1:5*nNodes);
% phi_dot = DV(5*nNodes+1:6*nNodes);
F       = DV(4*nNodes+1:5*nNodes);
theta   = DV(5*nNodes+1:6*nNodes);
m = DV(6*nNodes+1:7*nNodes);
tf      = DV(end);

dt = tf / N;

m0       = param.m;      % mass (kg) - change as needed
g       = param.g;     % lunar gravity (m/s^2)
I       = param.I;       % moment of inertia (kg*m^2)
d       = param.d;      % gimbal lever arm (m)
Fmax    = param.Fmax;  % max thrust (N)
Fmin    = param.Fmin;     % min thrust (N)
theta_max = param.theta_max; % rad
theta_min = param.theta_min;% rad


% --- compute accelerations at each node (vectorized) ---
% thrust direction in inertial frame = phi + theta
thrust_dir = theta;
xdd = (F./m) .* sin(thrust_dir);               % x double dot
zdd = (F./m) .* cos(thrust_dir) - g;           % z double dot
% phidd = (F .* d .* sin(theta)) ./ I;           % phi double dot (torque from offset)

% --- Collocation (trapezoidal) defects ---
% x dynamics: x_{k+1} - x_k - 0.5*dt*(x_dot_k + x_dot_{k+1}) = 0
ceq_x = x(2:end) - x(1:end-1) - 0.5*dt*( x_dot(1:end-1) + x_dot(2:end) );

% x_dot dynamics: x_dot_{k+1} - x_dot_k - 0.5*dt*(xdd_k + xdd_{k+1}) = 0
ceq_xdot = x_dot(2:end) - x_dot(1:end-1) - 0.5*dt*( xdd(1:end-1) + xdd(2:end) );

% z dynamics:
ceq_z = z(2:end) - z(1:end-1) - 0.5*dt*( z_dot(1:end-1) + z_dot(2:end) );

% z_dot dynamics:
ceq_zdot = z_dot(2:end) - z_dot(1:end-1) - 0.5*dt*( zdd(1:end-1) + zdd(2:end) );

% phi dynamics:
% ceq_phi = phi(2:end) - phi(1:end-1) - 0.5*dt*( phi_dot(1:end-1) + phi_dot(2:end) );
% 
% % phi_dot dynamics:
% ceq_phidot = phi_dot(2:end) - phi_dot(1:end-1) - 0.5*dt*( phidd(1:end-1) + phidd(2:end) );

m_dot = -F/(param.isp*9.81);
ceq_m = m(2:end) - m(1:end-1) - 0.5*dt*(m_dot(1:end-1) + m_dot(2:end));

x0       = BC.x0;
xdot0    = BC.xdot0;
z0       = BC.z0;
zdot0    = BC.zdot0;
% phi0     = BC.phi0;
% phidot0  = BC.phidot0;

xf       = BC.xf;
xdotf    = BC.xdotf;
zf       = BC.zf;
zdotf    = BC.zdotf;
% phif     = BC.phif;
% phidotf  = BC.phidotf;

ceq_bc = [ x(1)-x0;
           x_dot(1)-xdot0;
           z(1)-z0;
           z_dot(1)-zdot0;
           x(end)-xf;
           x_dot(end)-xdotf;
           z(end)-zf;
           z_dot(end)-zdotf
           m(1) - m0];

% --- collect all equality constraints ---
ceq = [ ceq_x(:); ceq_xdot(:); ceq_z(:); ceq_zdot(:); ceq_m(:); ceq_bc(:) ];

% thrust bounds (vector inequalities for each node)
c_F = [ F - Fmax;             % F <= Fmax
        Fmin - F ];           % F >= Fmin -> Fmin - F <= 0

% angle bounds
c_theta = [ theta - theta_max;
            theta_min - theta ];

% altitude non-negativity (z >= 0) -> -z <= 0
c_alt = -z+groundAlt(x); % Can add path contstraints here
% c_alt = -z;
% c_alt = 0;
c_mass = -m;
% stack all inequality constraints
c = [c_F(:); c_theta(:); c_alt(:); c_mass(:)];
end

function h = groundAlt(x)
n = length(x);
h = -10000*sin(x/5000-1.7)-10000;
end