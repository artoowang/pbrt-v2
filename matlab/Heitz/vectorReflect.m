% wi, wo: 3xN normalized vectors
% wh: 3x1 normalized vector
function wo = vectorReflect(wi, wh)
    N = size(wi, 2);
    wi_dot_wh = wh' * wi;  % 1xN
    if (sum(wi_dot_wh < 0) > 0)
        fprintf('There are wi(s) below wh''s hemisphere\n');
    end
    wo = 2 * repmat(wh, [1 N]) .* repmat(wi_dot_wh, [3 1]) - wi;
    wo(:, wi_dot_wh < 0) = 0;
end