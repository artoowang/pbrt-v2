% v: 3xN matrix
% threadDir: 3x1 vector
% threadN: 3x1 vector (assumed to be perpendicular to threadDir)
% theta, phi: 1xN vectors
function [theta, phi] = vector2threadLocal(v, threadDir, threadN)
    dotProductsT = min(max(threadDir' * v, -1), 1);   % 1xN
    theta = min(max(0.5*pi - acos(dotProductsT), -0.5*pi), 0.5*pi); % [-pi/2, pi/2]
    
    N = size(v, 2);
    vProj = v - repmat(dotProductsT, [3 1]) .* repmat(threadDir, [1 N]);
    dotProductsN = min(max(threadN' * vProj, -1), 1);   % 1xN
    
    [azi, elev] = cart2sph(v(1, :), v(2, :), v(3, :));
    theta = pi/2 - elev;
    phi = mod(azi, 2*pi);   % [0, 2*pi)
end