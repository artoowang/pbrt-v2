% w is in tangent space (i.e., geometric normal wg is (0,0,1))
% w is 3xM (M: number of vectors)
% density is 1xM
% Interates to one over projected hemisphere
% Note the density in the lower hemisphere is always 0
function density = D_Beckmann(w, alpha)
    cos_theta = w(3, :); % 1xM
    theta = acos(cos_theta); % 1xM
    density = exp(-(tan(theta)/alpha).^2) ./ (pi * alpha.^2 .* cos_theta.^4);
    density(cos_theta <= 0) = 0;
end