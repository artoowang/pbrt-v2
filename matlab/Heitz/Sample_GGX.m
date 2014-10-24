% wo: 3x1 vector
% wi: 3xN vectors
% pdf: 1xN vector
function [wi, pdf] = Sample_GGX(wo, N, alpha)
    u1 = rand(1, N);
    u2 = rand(1, N);
    theta_h = atan(alpha * sqrt(u1) ./ sqrt(1-u1));
    phi_h = u2 * 2 * pi;
    wh = sph2vector(theta_h, phi_h);
    cos_theta_h = wh(3, :);
    wi = vectorReflect(wo, wh);
    wo_dot_wh = wo' * wh;
    pdf = D_GGX(wh, alpha) .* cos_theta_h ./ (4 * wo_dot_wh);
    pdf(wo_dot_wh < 0) = 0;
end