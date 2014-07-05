% Currently only support single vector wi and wr (3x1 vector)
% Assume wi and wr in local tangent frame
function bsdfVal = AshikhminBSDF(wi, wr, densityFunction, f0)
    h = wi + wr;
    hLen = norm(h);
    if (hLen < 1e-6)
        bsdfVal = 0;
        return;
    end
    h = h / hLen;
    k_dot_h = wi' * h;
    % TODO: use gUsingGrid
    bsdfVal = densityFunction(h) * AshikhminAverageNH(densityFunction) ...
        * Fresnel(k_dot_h, f0) ...
        / (4 * Ashikhmin_g(wi, densityFunction) * Ashikhmin_g(wr, densityFunction));
end