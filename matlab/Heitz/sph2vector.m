% theta and phi will be columnized
% The result v has size 3xN, which N is the number of elements in theta and
% phi
function v = sph2vector(theta, phi)
    [x, y, z] = sph2cart(phi(:), pi/2-theta(:), 1);
    v = [x'; y'; z'];
end