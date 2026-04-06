% Read the image
Z = imread('Site04_final_adj_5mpp_surf.tif');

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

heightMap = Z;
% Assume 'heightMap' is your double array from the TIFF
minVal = min(heightMap(:));
maxVal = max(heightMap(:));

% 1. Shift the data so the lowest point is 0
normalizedMap = heightMap - minVal; 

% 2. Scale the data to fit the 16-bit range (0 to 65535)
% Check for divide-by-zero if the map is perfectly flat
range = maxVal - minVal;
if range > 0
    scaledMap = (normalizedMap / range) * 65535;
else
    scaledMap = normalizedMap;
end

% 3. Cast to uint16
finalMap = uint16(scaledMap);
imwrite(finalMap, 'LunarHeightmap_16bit.png', 'BitDepth', 16);
