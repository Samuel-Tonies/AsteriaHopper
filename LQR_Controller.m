m_l = 803; %kg
m_e = 150;
g = 1.625;
l = 4.2926/2;

% Based on linearized EOM about stablity point
A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, m_l*g/m_e, 0, 0;
     0, g*(m_e + m_l)/(l*m_e), 0, 0];

B = [0; 0; -l/m_e; 1/m_e];

% Arbitrartily
Q = [0.05, 0, 0, 0;
     0, 10000, 0, 0
     0, 0, .1, 0
     0, 0, 0, 50000];
R = 0.01;

K = lqr(A,B,Q,R);
