% wo: 3x1 vector
% wi: 3xN vectors
% pdf: 1xN vector
function pdf = Pdf_GGX(wo, wi, alpha)
    wh = halfVector(wi, wo);
    cos_theta_h = wh(3, :);
    wo_dot_wh = wo' * wh;
    if (any(wo_dot_wh < 0))
        error('dot(wo, wh) < 0?');
    end
    pdf = D_GGX(wh, alpha) .* cos_theta_h ./ (4 * wo_dot_wh);
    if (any(pdf < 0))
        error('pdf < 0?');
    end
    %pdf(wo_dot_wh < 0) = 0; % wo_dot_wh should be always >= 0
end