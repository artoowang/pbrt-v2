function vals = f_diffuse(wis)
    vals = ones(1, size(wis, 2)) / pi;
    vals(wis(3, :) < 0) = 0;
end