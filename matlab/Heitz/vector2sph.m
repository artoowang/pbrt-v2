% v: 3xN matrix
% theta, phi: 1xN vectors
function [theta, phi] = vector2sph(v)
    [azi, elev] = cart2sph(v(1, :), v(2, :), v(3, :));
    theta = pi/2 - elev;
    phi = mod(azi, 2*pi);   % [0, 2*pi)
end