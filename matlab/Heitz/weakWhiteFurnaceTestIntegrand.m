% wo: 3x1 vector
% wi: 3xN vectors
function val = weakWhiteFurnaceTestIntegrand(wi, wo, D, G1)
    wg_dot_wo = wo(3);
    if (wg_dot_wo < 0)
        error('wo shouldn''t be below horizon');
    end

    wh = halfVector(wi, wo); % 3xN
    densities = D(wh); % 1xN
    G1_vals = G1(wo, wh, D); % 1xN

    val = densities .* G1_vals ./ (4*wg_dot_wo);
end