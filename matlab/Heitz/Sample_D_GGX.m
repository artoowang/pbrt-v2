% w: 3xN vectors
% pdf: 1xN vector
function [w, pdf] = Sample_D_GGX(N, alpha)
    u1 = rand(1, N);
    u2 = rand(1, N);
    theta = atan(alpha * sqrt(u1) ./ sqrt(1-u1));
    phi = u2 * 2 * pi;
    w = sph2vector(theta, phi);
    cos_theta = w(3, :);
    pdf = D_GGX(w, alpha) .* cos_theta;
end