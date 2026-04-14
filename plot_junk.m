data = out;
X_i = data.simout.X_i;

z = X_i.data(:,3);
x = X_i.data(:,1);
t = X_i.Time;
g
euler_a_rad = data.simout.euler_a_rad.Data;
pitch = rad2deg(euler_a_rad(:,2));

figure
plot(t,z)
title("Vehicle Altitude")
xlabel("Time [s]")
ylabel("Altitude [m]")
grid on

figure
plot(t,pitch)
title("Vehicle Pitch")
xlabel("Time [s]")
ylabel("Pitch [deg]")
grid on