function trajPlot(sol,N,param)
%TRAJPLOT Plot results from Trajectory Tool

% Unpack results
x = sol(1:N+1);
xdot = sol(N+2:2*(N+1));
z = sol(2*(N+1)+1:3*(N+1));
zdot = sol(3*(N+1)+1:4*(N+1));
F = sol(4*(N+1)+1:5*(N+1));
theta = sol(5*(N+1)+1:6*(N+1));
mass_kg = sol(6*(N+1)+1:7*(N+1));
tf = sol(end);
t = linspace(0,tf,N+1);

figure
plot(t,x)
hold on
plot(t,z)
legend("X","Z")
grid on
title("X,Z vs Time")
xlabel("Time [s]")
ylabel("Distance X,Z [m]")

figure
plot(t,xdot)
hold on
plot(t,zdot)
legend("$\dot{X}$","$\dot{Z}$",'Interpreter','latex')
grid on
title("Velocity in X,Z vs Time")
xlabel("Time [s]")
ylabel("Velocity [m/s]")

figure
plot(t,F)
title("Thrust vs Time")
grid on
xlabel("Time [s]")
ylabel("Thrust [N]")

figure
plot(t,mass_kg)
grid on
title("Total Vehicle Mass vs Time")
xlabel("Time [s]")
ylabel("Mass [kg]")

figure
plot(x,z)
grid on
title("Trajectory Profile")
xlabel("Horizontal Distance [m]")
ylabel("Vertical Distance [m]")
if param.ground
    hold on
    plot(x,ground_alt(x))
    legend("Flight Trajectory","Ground")
end

end

