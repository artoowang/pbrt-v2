% wo: 3x1 vector
% wis: 3xN vectors
% D: @(w)
% G1: @(wo, wn, D)
% F: @(costheta)
% vals: 1xN vector
function vals = f_G1(wis, wo, D, G1, F)
    N = size(wis, 2);
    vals = zeros(1, N);
    
    wg_dot_wo = wo(3);
    if (wg_dot_wo < 0)
        error('wo shouldn''t be below horizon');
    elseif (wg_dot_wo == 0)
        return
    end
    
    abs_wg_dot_wis = abs(wis(3, :));
    
    wh = halfVector(wis, wo); % 3xN
    D_vals = D(wh); % 1xN
    G1_vals = G1(wo, wh, D); % 1xN
    F_vals = F(wo' * wh); % 1xN
    
    vals = F_vals .* G1_vals .* D_vals ./ abs_wg_dot_wis / wg_dot_wo / 4;
    
    vals(abs_wg_dot_wis == 0) = 0;
end