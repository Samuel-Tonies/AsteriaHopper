bim = blockedImage('ldem_87s_5mpp.tif');
sz = bim.Size; 

% 1. Define your "Resolution"
% Every 200x200 pixel block will become 1 single pixel in Z_small
blockSize = [200, 200]; 

% 2. Run the apply function
% This creates a new, much smaller blockedImage on your disk/memory
newBim = apply(bim, @(b) mean(b.Data, 'all', 'omitnan'), ...
               'BlockSize', blockSize, ...
               'UseParallel', true,'OutputLocation','test10');

Z_small = gather(newBim);

Z_small = double(Z_small);
[rows, cols] = size(Z_small);
[X, Y] = meshgrid(1:cols, 1:rows);

% Scale axes back to the original pixel units
X = X * blockSize(2);
Y = Y * blockSize(1);

figure;
surf(X, Y, Z_small, 'EdgeColor', 'none');
view(3); colorbar; axis tight; camlight(0,70);