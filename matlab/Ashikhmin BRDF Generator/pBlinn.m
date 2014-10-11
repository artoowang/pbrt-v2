% h_local is in tangent space (i.e., shading normal N is (0,0,1))
% h_local is 3xM (M: number of half vectors)
% density is 1xM
% pBlinn interates to one over hemisphere (for large enough exponent)
function density = pBlinn(h_local, exponent)
    n_dot_h = h_local(3,:);
    density = (exponent+1) .* n_dot_h.^exponent ./ (2*pi);
end