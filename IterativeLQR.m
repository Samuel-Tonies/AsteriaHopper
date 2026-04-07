function K = IterativeLQR(State)
coder.extrinsic('lqr');
K = zeros(1,4);
% Paramters assumed constant. Make sure match model inputs
m_l = 803; %kg
m_e = 150;
g = 1.625;
l = 4.2926/2;
% 
% x = State(1);
% phi = State(2);
% x_dot = State(3);
% phi_dot = State(4);
x = 0;
phi = 0;
x_dot = 0;
phi_dot = 0;

A = [0,                                                                                                                                                                                                                                                                                                                                        0, 1,                                                         0;
     0,                                                                                                                                                                                                                                                                                                                                        0, 0,                                                         1;
     0,                                                                                                                                    - (g*m_l*sin(phi)^2 - g*m_l*cos(phi)^2 + l*m_l*phi_dot^2*cos(phi))/(m_l*sin(phi)^2 + m_e) - (2*m_l*cos(phi)*sin(phi)*(- l*m_l*sin(phi)*phi_dot^2 + g*m_l*cos(phi)*sin(phi)))/(m_l*sin(phi)^2 + m_e)^2, 0,        -(2*l*m_l*phi_dot*sin(phi))/(m_l*sin(phi)^2 + m_e);
     0, (g*cos(phi))/l - (cos(phi)*(g*m_l*sin(phi)^2 - g*m_l*cos(phi)^2 + l*m_l*phi_dot^2*cos(phi)))/(l*(m_l*sin(phi)^2 + m_e)) - (sin(phi)*(- l*m_l*sin(phi)*phi_dot^2 + g*m_l*cos(phi)*sin(phi)))/(l*(m_l*sin(phi)^2 + m_e)) - (2*m_l*cos(phi)^2*sin(phi)*(- l*m_l*sin(phi)*phi_dot^2 + g*m_l*cos(phi)*sin(phi)))/(l*(m_l*sin(phi)^2 + m_e)^2), 0, -(2*m_l*phi_dot*cos(phi)*sin(phi))/(m_l*sin(phi)^2 + m_e)];

B = [0;
     0;
     1/(m_l*sin(phi)^2 + m_e);
     cos(phi)/(l*(m_l*sin(phi)^2 + m_e))];

% Arbitrartily MAKE SURE UPDATEd
Q = [10, 0, 0, 0;
     0, 200, 0, 0
     0, 0, .1, 0
     0, 0, 0, 1];

R = 0.8;

K = lqr(A,B,Q,R);

end

