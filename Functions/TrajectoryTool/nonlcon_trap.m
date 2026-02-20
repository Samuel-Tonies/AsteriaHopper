function [c, ceq] = nonlcon_trap(DV,N,param,BC)
% NONLCON_TRAP  nonlinear constraints (trapezoidal collocation) for 2D hopper
nNodes = N + 1;
x       = DV(1:nNodes);
x_dot   = DV(nNodes+1:2*nNodes);
z       = DV(2*nNodes+1:3*nNodes);
z_dot   = DV(3*nNodes+1:4*nNodes);
F       = DV(4*nNodes+1:5*nNodes);
theta   = DV(5*nNodes+1:6*nNodes);
m = DV(6*nNodes+1:7*nNodes);
tf      = DV(end);

dt = tf / N;

m0       = param.m;      % mass (kg) - change as needed
g       = param.g;     % lunar gravity (m/s^2)
Fmax    = param.Fmax;  % max thrust (N)
Fmin    = param.Fmin;     % min thrust (N)
theta_max = param.theta_max; % rad
theta_min = param.theta_min;% rad


% Thrust direction
thrust_dir = theta;
xdd = (F./m) .* sin(thrust_dir);               % x double dot
zdd = (F./m) .* cos(thrust_dir) - g;           % z double dot

% --- Trapezoidal Collocation ---
ceq_x = x(2:end) - x(1:end-1) - 0.5*dt*( x_dot(1:end-1) + x_dot(2:end) );
ceq_xdot = x_dot(2:end) - x_dot(1:end-1) - 0.5*dt*( xdd(1:end-1) + xdd(2:end) );
ceq_z = z(2:end) - z(1:end-1) - 0.5*dt*( z_dot(1:end-1) + z_dot(2:end) );
ceq_zdot = z_dot(2:end) - z_dot(1:end-1) - 0.5*dt*( zdd(1:end-1) + zdd(2:end) );

m_dot = -F/(param.isp*9.81);
ceq_m = m(2:end) - m(1:end-1) - 0.5*dt*(m_dot(1:end-1) + m_dot(2:end));

% Unpack BC
x0       = BC.x0;
xdot0    = BC.xdot0;
z0       = BC.z0;
zdot0    = BC.zdot0;
xf       = BC.xf;
xdotf    = BC.xdotf;
zf       = BC.zf;
zdotf    = BC.zdotf;

% BC Equalities
ceq_bc = [ x(1)-x0;
           x_dot(1)-xdot0;
           z(1)-z0;
           z_dot(1)-zdot0;
           x(end)-xf;
           x_dot(end)-xdotf;
           z(end)-zf;
           z_dot(end)-zdotf
           m(1) - m0];

% --- Collect all equality constraints ---
ceq = [ ceq_x(:); ceq_xdot(:); ceq_z(:); ceq_zdot(:); ceq_m(:); ceq_bc(:) ];

% Thrust bounds (vector inequalities for each node)
c_F = [ F - Fmax;             % F <= Fmax
        Fmin - F ];           % F >= Fmin -> Fmin - F <= 0

% Angle bounds
c_theta = [ theta - theta_max;
            theta_min - theta ];

% Altitude non-negativity (z >= 0) -> -z <= 0
if param.ground
    c_alt = -z + ground_alt(x);% Can add path contstraints here
else
    c_alt = 0;
end
c_mass = -m;
% Stack all inequality constraints
c = [c_F(:); c_theta(:); c_alt(:); c_mass(:)];
end