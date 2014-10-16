% wo: 3x1 vector
% wn: 3xN vectors. Note wn is only used to mask out the facets facing away
%     from wo, and to determine the result size
% result: 1xN vector
function result = G1_Beckmann(wo, wn, alpha)
    cos_theta_o = wo(3);
    theta_o = acos(cos_theta_o);
    a = 1 / (alpha * tan(theta_o));
    if a < 1.6
        lambda = (1 - 1.259*a + 0.396*a^2) / (3.535*a + 2.181*a^2);
    else
        lambda = 0;
    end
    result = 1 / (1 + lambda);
    result = repmat(result, [1 size(wn, 2)]);
    wo_dot_wn = wo' * wn; % 1xN
    result(wo_dot_wn < 0) = 0;
end