% w is in tangent space (i.e., geometric normal wg is (0,0,1))
% w is 3xM (M: number of vectors)
% density is 1xM
% D_Blinn interates to one over projected hemisphere
% Note the density in the lower hemisphere is always 0
function density = D_Blinn(w, exponent)
    n_dot_h = w(3, :);
    density = (exponent+2) .* n_dot_h.^exponent ./ (2*pi);
    density(n_dot_h < 0) = 0;
end