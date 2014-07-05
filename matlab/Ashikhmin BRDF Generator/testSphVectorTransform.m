function testSphVectorTransform()
    thetas = linspace(0, pi, 1000);
    phis = linspace(0, 2*pi, 1000);
    [pp, tt] = meshgrid(phis, thetas);
    V = sph2vector(tt(:), pp(:));
    [tt2, pp2] = vector2sph(V);
    maxError = max([abs(tt2(:)-tt(:)); abs(pp2(:)-pp(:))]);
    fprintf('maxError: %e\n', maxError);
end