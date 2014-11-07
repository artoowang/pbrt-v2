function val = Fresnel(k_dot_n, f0)
    val = f0 + (1-f0).*(1-k_dot_n).^5;
    val(k_dot_n < 0) = 0;
end
