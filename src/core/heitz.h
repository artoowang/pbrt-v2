
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

#ifndef PBRT_CORE_HEITZ_H
#define PBRT_CORE_HEITZ_H

#include "reflection.h"

// All the wo, wi, wh are in local tangent frame - shading normal is (0, 0, 1)
// Note in PBRT's implementation, the wo passed to BxDF::Sample_f() and wo, wi 
// passed to BxDF::f() are not guaranteed to be at the same side of the shading normal,
// i.e., their z component may be negative. To make things simpler (and matching other
// renderer better), my Heitz (BxDF) implementation handles this in its f() and 
// Sample_f(), so all later functions can assume wo is always at the same side as the
// shading normal (i.e., z component is non-negative).
// Notice wi can still be in the opposite side as wo - this is because PBRT determines 
// whether we're reflecting or transmitting using geometric normal, not shading normal.

// Integrate[D(wh)*cos(theta_h), {wh in hemisphere}] = 1
// D(wh) should be gauranteed to be zero in the lower hemisphere
// Integrate[Pdf(wo,wi), {wi in sphere}] = 1
// Geometric factor G is implemented here because it depends on D
class HeitzDistribution
{
public:
    HeitzDistribution() { }
    virtual ~HeitzDistribution() { }
    virtual float D(const Vector &wh) const = 0;
    virtual void Sample_f(const Vector &wo, Vector *wi,
                          float u1, float u2, float *pdf) const = 0;
    virtual float Pdf(const Vector &wo, const Vector &wi) const = 0;
    virtual float G(const Vector &wo, const Vector &wi, 
                    const Vector &wh) const = 0;
};

class GGXForHeitz : public HeitzDistribution
{
public:
    GGXForHeitz(float alpha);
    virtual float D(const Vector &wh) const;
    virtual void Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wo, const Vector &wi) const;
    virtual float G(const Vector &wo, const Vector &wi, const Vector &wh) const;

private:
    float mAlpha;
};

class Heitz : public BxDF
{
public:
    Heitz(const Spectrum &reflectance, Fresnel *f,
          const HeitzDistribution &d);
    Spectrum f(const Vector &woInput, const Vector &wiInput) const;
    Spectrum Sample_f(const Vector &wo, Vector *wi,
                      float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const;

    // TODO: for test
    bool mUseUniformSampling;

private:
    Spectrum R;
    const HeitzDistribution &mDistribution;
    Fresnel *fresnel;
};

#endif // PBRT_CORE_HEITZ_H
