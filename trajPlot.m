function trajPlot(sol,N,param)
clf;
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
% ylo = -4000;
% yhi = 500;
% patch([0 200 200 0],[ylo ylo yhi yhi], [.83 .51 .39],'FaceAlpha',.15,'EdgeColor','none')
% patch([200 3960 3960 200],[ylo ylo yhi yhi], [.21 .54 .87],'FaceAlpha',.15,'EdgeColor','none')
% patch([3960 7850 7850 3960],[ylo ylo yhi yhi], [.39 .60 .13],'FaceAlpha',.15,'EdgeColor','none')
% patch([7850 8000 8000 7850],[ylo ylo yhi yhi], [0.80 0.20 0.60],'FaceAlpha',.15,'EdgeColor','none')
% xlim([0 8000])
% text(5,   300, {'Phase 1', 'Launch'},     'HorizontalAlignment', 'center', 'FontSize', 10, 'Color','r')
% text(2000,  300, {'Phase 2', 'Coast'},      'HorizontalAlignment', 'center', 'FontSize', 10, 'Color','r')
% text(6000, 300, {'Phase 3', 'Mid-course'}, 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color','r')
% text(7900, 300, {'Phase 4', 'Landing'},    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color','r')
% annotation('textbox', [0.09 0.92 0.08 0.05], 'String', {'Phase 1','Ascent'}, ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10)
% annotation('textbox', [0.31 0.92 0.08 0.05], 'String', {'Phase 2','Ballistic Hop'}, ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10)
% annotation('textbox', [0.71 0.92 0.08 0.05], 'String', {'Phase 3','Pitch back'}, ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10)
% annotation('textbox', [0.85 0.92 0.08 0.05], 'String', {'Phase 4','Powered Descent'}, ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10)

figure
plot(t,rad2deg(theta))
grid on
title("Pitch vs Time")
xlabel("Time [s]")
ylabel("Pitch [deg]")

end

