% k is assumed to be in the tangent frame (where surface normal n is +Z)
% Currently we only support a single 3x1 vector k
function val = Ashikhmin_g(k, densityFunction)
    % Find the transform from vector n to vector k: it's a rotation along 
    % the vector cross(n, k). This rotation transform a vector in the
    % hemisphere defined by k, to the hemishpere defined by n
    rotation_vec = vrrotvec([0 0 1]', k);
    R = vrrotvec2mat(rotation_vec);

    val = quad2d(@(theta, phi) integrand(theta, phi, R, densityFunction), ...
        0, pi/2, 0, 2*pi);
end

function val = integrand(theta, phi, R, densityFunction)
    % Note: the local space (theat, phi) here are defined by k, 
    %       not n (surface normal). The transform from k space to n space
    %       is given by R
    h = sph2vector(theta, phi); % in k space, 3xN, (N: number of thetas)
    % By definition, dot(h,k) = cos(theta) = h(3)
    h_dot_k = h(3, :);  % 1xN
    h_n = R * h;    % in n space, 3xN
    densities = densityFunction(h_n);   % 1xN
    densities = reshape(densities(:), size(theta));
    h_dot_k = reshape(h_dot_k(:), size(theta));
    % N dot H is cos(theta) by definition
    val = h_dot_k .* densities .* sin(theta);
end