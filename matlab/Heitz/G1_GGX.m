% wo: 3x1 vector
% wn: 3xN vectors. Note wn is only used to mask out the facets facing away
%     from wo, and to determine the result size
% result: 1xN vector
function result = G1_GGX(wo, wn, alpha)
    cos_theta_o = wo(3);
    theta_o = acos(cos_theta_o);
    % When theta_o == 0, a == Inf, and later 1/a^2 == 0
    a = 1 / (alpha * tan(theta_o));
    lambda = (-1 + sqrt(1 + 1 / a^2)) / 2;
    result = 1 / (1 + lambda);
    result = repmat(result, [1 size(wn, 2)]);
    wo_dot_wn = wo' * wn; % 1xN
    result(wo_dot_wn < 0) = 0;
end