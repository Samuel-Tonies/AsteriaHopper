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

F_init = zeros(N+1,1);
for i=1:ceil((N+1)/2)
        F_init(i) = param.Fmax*(1-2*t_guess(i)/tf0);
end
for i=ceil((N+1)/2)+1:N+1
        F_init(i) = F_init(i-ceil((N+1)/2));
end
    
theta_init = atan2(zd_init,xd_init);
m_init = param.m - t_guess/tf0*prop_guess;

DV0 = [x_init;xd_init; z_init; zd_init; F_init;...
    theta_init; m_init; tf0];
end