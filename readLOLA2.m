% Read the image
Z = imread('LM1_final_adj_5mpp_surf.tif');

% % Convert to double for math operations
Z = double(Z);

% Create a grid the same size as the image
[X, Y] = meshgrid(1:size(Z,2), 1:size(Z,1));

% Generate the 3D surface
figure;
surf(X, Y, Z, 'EdgeColor', 'none'); % 'none' removes the black grid lines

% Improve the visuals
% colormap terrain;       % Uses a natural color map
colorbar;               % Shows the height scale
view(3);                % Sets the camera to a 3D perspective
axis tight;             % Removes empty space around the map
camlight(0,70);               % Adds a light source to show shadows
% lighting phong;         % Makes the surface look smooth
