
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
#include "materials/heitz.h"
#include "spectrum.h"
#include "core/heitz.h"
#include "paramset.h"
#include "texture.h"

// PlasticMaterial Method Definitions
BSDF *HeitzMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                             const DifferentialGeometry &dgShading,
                             MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (mBumpMap)
        Bump(mBumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    // TODO: implement
    //Spectrum kd = mKd->Evaluate(dgs).Clamp();
    //if (!kd.IsBlack()) {
    //    BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
    //    bsdf->Add(diff);
    //}

    //float etat = mEtaT->Evaluate(dgs);
    Spectrum ks = mKs->Evaluate(dgs).Clamp();
    if (!ks.IsBlack()) {
        // TODO: test
        //Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.f, etat); // Note this is different from uber, which I think implements this wrong
        Fresnel *fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
        float rough = mRoughness->Evaluate(dgs);

        Heitz *spec = NULL;

        const HeitzDistribution &distribution = *(BSDF_ALLOC(arena, GGXForHeitz)(rough));
        spec = BSDF_ALLOC(arena, Heitz)(ks, fresnel, distribution);

        bsdf->Add(spec);
    }
    return bsdf;
}

HeitzMaterial *CreateHeitzMaterial(const Transform &xform,
        const TextureParams &mp)
{
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<float> > etat = mp.GetFloatTexture("index", 1.5f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    string ndfFilePath = mp.FindString("ndffilepath", "");

    HeitzMaterial *mtl = new HeitzMaterial(Kd, Ks, roughness, etat, bumpMap, ndfFilePath);

    return mtl;
}

