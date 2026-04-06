function h = ground_alt(x)
% Function that returns the altitude/height of ground relative to starting
% point for a given position

% Suppose to mimick shackleton crater
% Function should be differentiable or else solver issues can occur

n = length(x);
h = -4000*sin(x/250-.7);
% A1 = 1;
% sig1 = 1;
% A2 = 1;
% sig2 = 1;
% h = -4000*sin(2*pi*x/42000) - A1*tanh(x/sig1) + A2*(tanh(x-21000/sig2)...
%     - tanh(-21000/sig2));
% h = (x-21000/2).^2 - 4100;
end