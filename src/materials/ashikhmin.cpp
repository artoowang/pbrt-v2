
/*
    pbrt source code Copyright(c) 2014 Chun-Po Wang.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


#include "stdafx.h"
#include "materials/ashikhmin.h"
#include "spectrum.h"
#include "core/ashikhmin.h"
#include "paramset.h"
#include "texture.h"

#include "mipmap.h"

// PlasticMaterial Method Definitions
BSDF *AshikhminMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                                 const DifferentialGeometry &dgShading,
                                 MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum ks = Ks->Evaluate(dgs).Clamp();
    if (!ks.IsBlack()) {
        // TODO: test
        //Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
        Fresnel *fresnel = BSDF_ALLOC(arena, FresnelNoOp);
        float rough = roughness->Evaluate(dgs);
        BxDF *spec = BSDF_ALLOC(arena, Ashikhmin)
                       (ks, fresnel, BSDF_ALLOC(arena, BlinnForAshikhmin)(1.f / rough));
        bsdf->Add(spec);
    }
    return bsdf;
}

// s, t in [0, 1]
float
AshikhminMaterial::testMIPMapFunc(float s, float t)
{
    const float a = 2.f, b = 3.f, c = 4.f;
    return a * s + b * t + c;
}

void
AshikhminMaterial::testMIPMap(void)
{
    // Currently this only works when width and height are power of 2
    // (i.e., the maxErr is about 1e-7 ~ 1e-6)
    // Otherwise, the resampling method seems creates error (around 1e-2)
    const int width = 3, height = 3;
    float *map = new float[width * height];

    // The center of the upper-left pixel (x,y) = (0,0) is mapped to
    // (s,t) = (0,0)+half_pixel_size, while the lower-right pixel (x,y) = (width-1,height-1)
    // is mapped to (s,t) = (1,1)-half_pixel_size
    // I.e., s = (x+0.5)/width, t = (y+0.5)/height
    for (int y = 0; y < height; ++y) {
        const float t = (y + 0.5f) / height;
        for (int x = 0; x < width; ++x) {
            const float s = (x + 0.5f) / width;
            map[y*width + x] = testMIPMapFunc(s, t);
        }
    }

    // doTrilinear and maxAnisotropy should make no difference if we only use basic Lookup()
    MIPMap<float> mipmap(width, height, map, false, 8.f, TEXTURE_CLAMP);

    // TODO: test print map and the resampled map
    // TODO: yes, the resampling indeed does some strange things. Need to check the PBRT book
    fprintf(stderr, "map:\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(stderr, " %.2f", map[y*width + x]);
        }
        fprintf(stderr, "\n");
    }
    const int testWidth2 = 9, testHeight2 = 9;
    fprintf(stderr, "resampled map:\n");
    for (int y = 0; y < testHeight2; ++y) {
        const float t = (y + 0.5f) / testHeight2;
        for (int x = 0; x < testWidth2; ++x) {
            const float s = (x + 0.5f) / testWidth2;
            fprintf(stderr, " %.2f", mipmap.Lookup(s, t));
        }
        fprintf(stderr, "\n");
    }

    // Test the area between [0.5/width,1-0.5/width]x[0.5/height,1-0.5/height]
    // (1-1/width by 1-1/height)
    float maxErr = 0.f;
    const int testWidth = 100, testHeight = 100;
    for (int y = 0; y < testHeight; ++y) {
        const float t = 0.5f/height + (1-1.f/height) * y / (testHeight-1);
        for (int x = 0; x < testWidth; ++x) {
            const float s = 0.5f/width + (1-1.f/width) * x / (testWidth-1);
            const float funcVal = testMIPMapFunc(s, t),
                        mipmapVal = mipmap.Lookup(s, t),
                        err = fabsf(funcVal - mipmapVal);
            maxErr = max(maxErr, err);
        }
    }
    fprintf(stderr, "MIPMap test: max error: %e\n", maxErr);

    delete [] map;
}

AshikhminMaterial *CreateAshikhminMaterial(const Transform &xform,
        const TextureParams &mp)
{
    // TODO: tests
    //Ashikhmin::testSphVectorTransform();
    Ashikhmin::testAverageNHAndFactor_g();
    //Quaternion::unitTest();
    //AshikhminMaterial::testMIPMap();

    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");

    // TODO: test
    AshikhminCache::sGridResolution = mp.FindInt("gridresolution", 32);

    return new AshikhminMaterial(Ks, roughness, bumpMap);
}


