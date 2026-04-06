bim = blockedImage('ldem_87s_5mpp.tif');
sz = bim.Size; 

% 1. Define your "Resolution"
% Every 200x200 pixel block will become 1 single pixel in Z_small
blockSize = [120, 120]; 

% 2. Run the apply function
% This creates a new, much smaller blockedImage on your disk/memory
newBim = apply(bim, @(b) mean(b.Data, 'all', 'omitnan'), ...
               'BlockSize', blockSize, ...
               'UseParallel', true,'OutputLocation','lidarBlocked7');

Z_small = gather(newBim);

Z_small = double(Z_small);
[rows, cols] = size(Z_small);
[X, Y] = meshgrid(1:cols, 1:rows);

% Scale axes back to the original pixel units
X = X * blockSize(2);
Y = Y * blockSize(1);

figure;
surf(X, Y, Z_small, 'EdgeColor', 'none');
view(3); cd = colorbar; axis tight; camlight(0,70);
cd.Color = 'w';
axis off;
grid off;

% 2. Match the background color 
% Choose 'w' for white, 'k' for black, or a custom [R G B]
bgColor = 'k'; 
set(gcf, 'Color', bgColor); % Sets the window background
set(gca, 'Color', bgColor); % Sets the plot area background


% 1. Define Start and End Points (in pixel units)
% Adjust these values to place the arc where you want on your map
startX = 19458;  startY = 21804; 
endX = 20907;    endY = 21252;

% 2. Define the Peak Height (how high the arc goes above the ground)
archHeight = 1250 

% 3. Create the horizontal path
numPoints = 100; % Smoothness of the line
trajX = linspace(startX, endX, numPoints);
trajY = linspace(startY, endY, numPoints);

% 4. Calculate the Ground Elevation under the path
% We interpolate Z_small to ensure the arc starts/ends at the surface
startZ = interp2(X, Y, Z_small, startX, startY);
endZ = interp2(X, Y, Z_small, endX, endY);

% 5. Calculate the Parabola (Z-axis)
% A simple normalized parabola: y = -4*h*(x-0.5)^2 + h
t = linspace(0, 1, numPoints);
relativeHeight = -4 * archHeight * (t - 0.5).^2 + archHeight;
% Add the linear slope between start and end elevation
trajZ = linspace(startZ, endZ, numPoints) + relativeHeight;

% 6. Plot the trajectory
hold on;
% We'll use a thick, glowing line (Neon Cyan) to pop against the black background
plot3(trajX, trajY, trajZ, 'Color', [1 0 0], 'LineWidth', 3);

% Optional: Add a "landing target" or start point marker
plot3(startX, startY, startZ, 'wx', 'MarkerSize', 10);
plot3(endX, endY, endZ, 'wx', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
