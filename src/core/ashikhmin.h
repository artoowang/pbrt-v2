
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_ASHIKHMIN_H
#define PBRT_CORE_ASHIKHMIN_H

#include "reflection.h"

// D is normalized for hemisphere (instead of the projected hemisphere as in Blinn)
// I.e., Integrate[D(wh), {wh in hemisphere}] = 1
// Currently Sample_f() and Pdf() implementation is the same as Blinn
// I think Pdf() is against wi instead of wh, but still need to confirm that
class BlinnForAshikhmin : public MicrofacetDistribution {
public:
    BlinnForAshikhmin(float e);
    float D(const Vector &wh) const;
    virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wi, const Vector &wo) const;
private:
    float exponent;
};

class Ashikhmin : public BxDF {
public:
    Ashikhmin(const Spectrum &reflectance, Fresnel *f,
        MicrofacetDistribution *d);
    Spectrum f(const Vector &wo, const Vector &wi) const;
    Spectrum Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const;

    // For tests
    static void testSphVectorTransform(void);
    static void testAverageNHAndFactor_g(void);

private:
    float averageNH(void) const;
    float gFactor(const Vector &v) const;

    static int averageNHIntegrand(const int *ndim, const double xx[],
            const int *ncomp, double ff[], void *userdata);
    static int gFactorIntegrand(const int *ndim, const double xx[],
                const int *ncomp, double ff[], void *userdata);

    Spectrum R;
    MicrofacetDistribution *distribution;
    Fresnel *fresnel;

    struct gFactorIntegrandData {
        const MicrofacetDistribution *distribution;
        Transform nToV;
    };
};

#endif // PBRT_CORE_ASHIKHMIN_H
