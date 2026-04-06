data = out;
X_i = data.simout.X_i;

z = X_i.data(:,3);
x = X_i.data(:,1);
t = X_i.Time;

figure
plot(t,x)
hold on
plot(t,z)
title("Vehicle States")