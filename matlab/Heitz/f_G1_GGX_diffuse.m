function vals = f_G1_GGX_diffuse(wis, wo, alpha, F)
    vals = f_G1_diffuse(wis, wo, @(w) D_GGX(w, alpha), @(wo, wn, D) G1_GGX(wo, wn, alpha), F);
end