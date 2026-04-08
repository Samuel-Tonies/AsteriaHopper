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
Q = [.00001, 0, 0, 0;
     0, 100, 0, 0
     0, 0, .1, 0
     0, 0, 0, 1];
R = 0.8;

K = lqr(A,B,Q,R);
Acl = A - B*K;
eig(Acl)

K_x = K;
K_y = [K(1), -K(2), K(3), -K(4)];
