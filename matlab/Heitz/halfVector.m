% wi: 3xN normalized vectors
% wo: 3x1 normalized vector
function wh = halfVector(wi, wo)
    N = size(wi, 2);
    wh = wi + repmat(wo, [1 N]);
    lens = sqrt(sum(wh .* wh, 1));
    wh = wh ./ repmat(lens, [3 1]);
end