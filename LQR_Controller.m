m_l = 803; %kg
m_e = 150;
g = 1.625;
l = 4.2926/2;

% Based on linearized EOM about stablity point
A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, m_l*g/m_e, 0, 0;
     0, g*(m_e + m_l)/(l*m_e), 0, 0];

B = [0; 0; 1/m_e; 1/(m_e*l)];

% Arbitrartily
Q = [0.1, 0, 0, 0;
     0, 200, 0, 0
     0, 0, .1, 0
     0, 0, 0, 1];
R = .1;

K = lqr(A,B,Q,R);
Acl = A - B*K;
eig(Acl)

% syms x phi x_dot phi_dot T real % Define state variables
% x_ddot = (m_l*g*sin(phi)*cos(phi) - m_l*l*phi_dot^2*sin(phi)) / (m_e+m_l*(sin(phi)^2));
% phi_ddot = x_ddot*cos(phi)/l + g/l*sin(phi);
% 
% A_nonl = [x_dot;
%      phi_dot;
%      x_ddot;
%      phi_ddot];
% 
% A_jacob = jacobian(A_nonl,[x;phi;x_dot;phi_dot]);
% 
% B_nonl = [0;
%           0;
%           T/(m_e+m_l*(sin(phi)^2));
%           T*cos(phi)/(l*(m_e+m_l*(sin(phi)^2)))];
% B_jacob = jacobian(B_nonl, T);