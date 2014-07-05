% Currently only support single 3x1 vector k
% k is in local tangent space (where n is (0, 0, 1))
function P = AshikhminP(k, densityFunction)
    n_dot_k = k(3);
    P = n_dot_k * AshikhminAverageNH(densityFunction) ...
        / Ashikhmin_g(k, densityFunction);
end