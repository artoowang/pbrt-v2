function saveNDF(NDF, filepath)
    fp = fopen(filepath, 'w');
    [thetaRes, phiRes] = size(NDF);
    fprintf(fp, '%d %d\n', thetaRes, phiRes);
    for t = 1:thetaRes
        for p = 1:phiRes
            fprintf(fp, '%e ', NDF(t, p));
        end
        fprintf(fp, '\n');
    end
    fclose(fp);
end