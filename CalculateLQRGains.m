% function K = CalculateLQRGains(States)  
% Paramters assumed constant. Make sure match model inputs
% m_l = 803; %kg
% m_e = 150;
% g = 1.625;
% l = 4.2926/2;

% assume current state value is linearization point
% x0 = States(1);
% phi0 = States(2);
% x_dot0 = States(3);
% phi_dot0 = States(4);

syms x phi x_dot phi_dot T m_e m_l g l real % Define state variables
x_ddot = (m_l*g*sin(phi)*cos(phi) - m_l*l*phi_dot^2*sin(phi)) / (m_e+m_l*(sin(phi)^2));
phi_ddot = x_ddot*cos(phi)/l + g/l*sin(phi);

A_nonl = [x_dot;
     phi_dot;
     x_ddot;
     phi_ddot];

A_jacob = jacobian(A_nonl,[x;phi;x_dot;phi_dot]);
% A_lin = double(subs(A_jacob, [x;phi;x_dot;phi_dot], [x0;phi0;x_dot0;phi_dot0]));

B_nonl = [0;
          0;
          T/(m_e+m_l*(sin(phi)^2));
          T*cos(phi)/(l*(m_e+m_l*(sin(phi)^2)))];
B_jacob = jacobian(B_nonl, T);
B_lin = double(subs(B_jacob, [x;phi;x_dot;phi_dot], [x0;phi0;x_dot0;phi_dot0]));

% Arbitrartily
Q = [10, 0, 0, 0;
     0, 200, 0, 0
     0, 0, .1, 0
     0, 0, 0, 1];
R = 0.8;

K = lqr(A_lin,B_lin,Q,R);

% end

