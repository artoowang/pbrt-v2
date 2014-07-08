% Support single wr but vectorized wi
% (wr: 3x1 vector, wi: 3xN matrix)
% Assume wi and wr in local tangent frame
function bsdfVal = AshikhminBSDFUsingGrid(wi, wr, f0, gGrid, avgNH)
    h = [wi(1,:) + wr(1); wi(2,:) + wr(2); wi(3,:) + wr(3)];
    hLen = sqrt(sum(h.^2, 1));
    validInputs = hLen > 1e-6;
    
    bsdfVal = zeros(1, numel(validInputs));
    if (sum(validInputs) == 0)
        % No valid input
        return;
    end
    
    h = h(:, validInputs) ./ repmat(hLen(validInputs), [3 1]);
    k_dot_h = wr' * h;
    
    denominator = 4 .* Ashikhmin_gUsingGrid(wi, gGrid) ...
        .* Ashikhmin_gUsingGrid(wr, gGrid);
    bsdfValValid = gGrid.densityFunction(h) .* avgNH ...
        .* Fresnel(k_dot_h, f0) ./ denominator;
    bsdfVal(validInputs) = bsdfValValid;
end