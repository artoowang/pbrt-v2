
/*
    pbrt source code Copyright(c) 1998-2014 Chun-Po Wang.

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
#include "heitz.h"

#include <boost/format.hpp>

const float sSmallValue = 1e-6f;

// ----------------------------------------------------------------------------

GGXForHeitz::GGXForHeitz(float alpha) :
    mAlpha(alpha)
{
    // TODO: this might still be too small since we're using mAlpha^2
    if (mAlpha < sSmallValue) {
        mAlpha = sSmallValue;
    }
}

float
GGXForHeitz::D(const Vector &wh) const
{
    float costhetah = CosTheta(wh);
    if (costhetah < sSmallValue) {
        return 0.f;
    } else {
        Assert(costhetah > 0.f && costhetah <= 1.f);
        float thetah = acosf(costhetah),
              costhetah2 = costhetah * costhetah,
              costhetah4 = costhetah2 * costhetah2,
              tanthetah = tanf(thetah),
              alphaSq = mAlpha * mAlpha,
              tanOverAlpha2 = (tanthetah * tanthetah) / alphaSq,
              factor = 1.f + tanOverAlpha2;
        return 1.f / (M_PI * alphaSq * costhetah4 * factor * factor);
    }
}

// wo.z >= 0 is guaranteed
// Note the sampled wi ~ D(wh)*cos(thetah)/(4*dot(wo,wh))
// and pdf = 0 when dot(wo, wh) < 0 (make sure Pdf() does the same thing)
void
GGXForHeitz::Sample_f(const Vector &wo, Vector *wi, float u1, float u2,
                     float *pdf) const
{
    // We allow cos(theta_o) == 0 because it doesn't break our computation
    Assert(CosTheta(wo) >= 0.f);

    // From [Walter et al., 2007]
    float thetah = atanf(mAlpha * sqrt(u1) / sqrt(1.f-u1)),
          costhetah = cosf(thetah), 
          sinthetah = sinf(thetah),
          phih = u2 * 2.f * M_PI;
    Vector wh = SphericalDirection(sinthetah, costhetah, phih);

    // Compute incident direction by reflecting about $\wh$
    float dotHO = Dot(wo, wh);

    if (dotHO < sSmallValue) {
        *pdf = 0.f;

    } else {
        *wi = -wo + 2.f * dotHO * wh;
        // TODO: test
        //if (wi->z < 0) {
        //    fprintf(stderr, "wi.z == %f\n", wi->z);
        //}
        // TODO: can be optimized
        *pdf = D(wh) * costhetah / (4.f * dotHO);
    }
}

float
GGXForHeitz::Pdf(const Vector &wo, const Vector &wi) const
{
    // We allow cos(theta_o) == 0 because it doesn't break our computation
    Assert(CosTheta(wo) >= 0.f);

    Vector wh = Normalize(wo + wi);
    float dotHO = Dot(wo, wh);
    if (dotHO < sSmallValue) {
        return 0.f;
    }

    return D(wh) * CosTheta(wh) / (4.f * dotHO);
}

float
GGXForHeitz::G(const Vector &wo, const Vector &/*wi*/, const Vector &wh) const
{
    float costhetao = CosTheta(wo);
    if (costhetao < sSmallValue) {
        return 0.f;
    }

    // TODO: currently, we only compute G1(wo, wh)
    float dotHO = Dot(wo, wh);
    if (dotHO < 0.f) {
        return 0.f;
    }
    float thetao = acosf(costhetao);
    if (thetao < sSmallValue) {
        return 1.f;
    }
    float a = 1.f / (mAlpha * tanf(thetao)),
          lambda = (-1.f + sqrtf(1.f + 1.f / (a*a))) / 2;
    return 1.f / (1.f + lambda);
}

// ----------------------------------------------------------------------------

Heitz::Heitz(const Spectrum &reflectance, Fresnel *f,
             const HeitzDistribution &d)
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
      mUseUniformSampling(false),   // TODO: test
      R(reflectance), mDistribution(d), fresnel(f)
{
}

Spectrum
Heitz::f(const Vector &woInput, const Vector &wiInput) const
{
    // In PBRT's implementation, woInput and wiInput is not guaranteed to be at the same side of the shading normal
    // We inverse them if woInput is on the other side.
    // Note: because BSDF::f() determines BSDF_REFLECTION and BSDF_TRANSMISSION using ng (geometric normal) instead
	//       of nn (shading normal), it's possible that woInput and wiInput has opposite sign of z-axis. Therefore,
	//       we determine if the shading normal is at the opposite side using woInput (viewing direction) only: we
	//       always make the viewing direction from the +Z side

    Vector wo = woInput, wi = wiInput;
    if (wo.z < 0.f) {
        // Reverse side
        wo = -wo;
        wi = -wi;
    }

    float cosThetaO = CosTheta(wo),
          absCosThetaI = AbsCosTheta(wi);
    Assert(cosThetaO >= 0.f);
    if (cosThetaO < sSmallValue || absCosThetaI < sSmallValue) {
        return Spectrum(0.f);
    }

    Vector wh = wi + wo;
    if (wh.x == 0. && wh.y == 0. && wh.z == 0.) {
        return Spectrum(0.f);
    }
    wh = Normalize(wh);
    float i_dot_h = Dot(wi, wh);
    Spectrum F = fresnel->Evaluate(i_dot_h);

    // Specular/glossy component
    float spec = mDistribution.D(wh) * mDistribution.G(wo, wi, wh) /
            (4.f * cosThetaO * absCosThetaI);

    // Diffuse component
    float diff = INV_PI;

    // TODO: test
    //F = Spectrum(0.f);

    // Note wi is possible to be at lower hemisphere, so we need to use
    // AbsCosTheta()
    return R * (F * spec + (Spectrum(1.f) - F) * diff);

    // TODO: test
    /*float Dval = mDistribution.D(wh),
          Gval = mDistribution.G(wo, wi, wh);
    Spectrum val = R * Dval * Gval * F /
            (4.f * cosThetaO * absCosThetaI);
    if (isinf(val.y()) || val.y() < -sSmallValue) {
        fprintf(stderr, "invalid value\n");
    }
    return val;*/
}

Spectrum
Heitz::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const 
{
    // TODO: for test
    if (mUseUniformSampling) {
        return BxDF::Sample_f(wo, wi, u1, u2, pdf);

    } else {
        // TODO: what percentage between diffuse and specular? Currently use 50%
        float pdfSpec = 0.f, pdfDiff = 0.f;
        if (u1 < 0.5f) {
            // Stretch u1 back to [0, 1)
            u1 *= 2;

            if (wo.z < 0.f) {
                // Reverse side
                mDistribution.Sample_f(-wo, wi, u1, u2, &pdfSpec);
                *wi = -(*wi);
            } else {
                mDistribution.Sample_f(wo, wi, u1, u2, &pdfSpec);
            }

            // Note: wo and wi are not guaranteed to be on the +Z side, nor are
            //      they to be on the same side

            if (pdfSpec > 0.f) {
                // The sample is valid. Compute the PDF from diffuse sampling
                // method
                pdfDiff = SameHemisphere(wo, *wi) ? AbsCosTheta(*wi) * INV_PI : 0.f;

            } else {
                // The sample is invalid - note it still needs to be considered
                // in Monte Carlo because invalid sample takes part of the PDF
                *pdf = 0.f;
                return Spectrum(0.f);
            }

        } else {
            // Stretch u1 back to [0, 1)
            u1 = 2 * (u1 - 0.5f);

            // Sampling wi at the same side as wo
            *wi = CosineSampleHemisphere(u1, u2);
            if (wo.z < 0.f) {
                *wi = -(*wi);
            }

            // Compute PDF for diffuse
            Assert(SameHemisphere(wo, *wi));
            pdfDiff = AbsCosTheta(*wi) * INV_PI;

            // Compute PDF for specular
            if (wo.z < 0.f) {
                // Reverse side
                pdfSpec = mDistribution.Pdf(-wo, -(*wi));
            } else {
                pdfSpec = mDistribution.Pdf(wo, *wi);
            }
        }

        Assert(pdf);
        *pdf = (pdfSpec + pdfDiff) * 0.5f;
        return f(wo, *wi);
    }
}

float
Heitz::Pdf(const Vector &wo, const Vector &wi) const 
{
    // TODO: for test
    if (mUseUniformSampling) {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0.f;

    } else {
        // TODO: what percentage between diffuse and specular? Currently use 50%
        float pdfSpec = 0.f;
        if (wo.z < 0.f) {
            // Reverse side
            pdfSpec = mDistribution.Pdf(-wo, -wi);
        } else {
            pdfSpec = mDistribution.Pdf(wo, wi);
        }

        float pdfDiff = SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0.f;
        return (pdfDiff + pdfSpec) * 0.5f;
    }
}
