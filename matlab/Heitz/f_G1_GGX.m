function vals = f_G1_GGX(wis, wo, alpha, f0)
    vals = f_G1(wis, wo, @(w) D_GGX(w, alpha), @(wo, wn, D) G1_GGX(wo, wn, alpha), @(costheta) Fresnel(costheta, f0));
end