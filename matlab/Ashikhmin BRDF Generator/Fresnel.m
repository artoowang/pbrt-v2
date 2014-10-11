function val = Fresnel(k_dot_n, f0)
    k_dot_n(k_dot_n < 0) = nan;
    val = f0 + (1-f0).*(1-k_dot_n).^5;
end