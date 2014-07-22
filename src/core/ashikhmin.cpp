
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
#include "ashikhmin.h"

#include <sstream>

// For numerical integration
#include "cuba.h"

using std::stringstream;

BlinnForAshikhmin::BlinnForAshikhmin(float e)
{
    if (e > 10000.f || isnan(e)) {
        e = 10000.f;
    }
    exponent = e;
}

float
BlinnForAshikhmin::D(const Vector &wh) const
{
    // Note: unlike class Blinn, we want to make sure there is no microfacet facing downward
    float costhetah = max(CosTheta(wh), 0.f);    // TODO: Interesting, this breaks PBRT. Why?
    return (exponent+1) * INV_TWOPI * powf(costhetah, exponent);
}

// Currently BlinnForAshikhmin::Sample_f() and Pdf() is the same as Blinn
void
BlinnForAshikhmin::Sample_f(const Vector &wo, Vector *wi, float u1, float u2,
                     float *pdf) const
{
    // Compute sampled half-angle vector $\wh$ for Blinn distribution
    float costheta = powf(u1, 1.f / (exponent+1));
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    float phi = u2 * 2.f * M_PI;
    Vector wh = SphericalDirection(sintheta, costheta, phi);
    if (!SameHemisphere(wo, wh)) wh = -wh;

    // Compute incident direction by reflecting about $\wh$
    *wi = -wo + 2.f * Dot(wo, wh) * wh;

    // Compute PDF for $\wi$ from Blinn distribution
    float blinn_pdf = ((exponent + 1.f) * powf(costheta, exponent)) /
                      (2.f * M_PI * 4.f * Dot(wo, wh));
    if (Dot(wo, wh) <= 0.f) blinn_pdf = 0.f;
    *pdf = blinn_pdf;
}

float
BlinnForAshikhmin::Pdf(const Vector &wo, const Vector &wi) const
{
    Vector wh = Normalize(wo + wi);
    float costheta = AbsCosTheta(wh);
    // Compute PDF for $\wi$ from Blinn distribution
    float blinn_pdf = ((exponent + 1.f) * powf(costheta, exponent)) /
                      (2.f * M_PI * 4.f * Dot(wo, wh));
    if (Dot(wo, wh) <= 0.f) blinn_pdf = 0.f;
    return blinn_pdf;
}

string
BlinnForAshikhmin::signature(void) const
{
    stringstream ss;
    ss.precision(5);
    ss << std::scientific;
    ss << "BlinnForAshikhmin:" << exponent;
    return ss.str();
}


AshikhminCache::AshikhminCacheMap AshikhminCache::sCache;

AshikhminCache::AshikhminCache() :
        mDistribution(NULL), mAvgNH(1.f)
{
}

AshikhminCache::AshikhminCache(const MicrofacetDistribution &distribution) :
        mDistribution(&distribution),
        mAvgNH(Ashikhmin::computeAverageNH(distribution))
{
    // TODO: Init g grid
}

float
AshikhminCache::gFactor(const Vector &v) const
{
    // TODO
    return 1.f;
}

float
AshikhminCache::averageNH(void) const
{
    return mAvgNH;
}

const AshikhminCache&
AshikhminCache::get(const MicrofacetDistribution &distribution)
{
    // TODO: need to use mutex to make it thread safe
    AshikhminCacheMap::const_iterator it;
    if ((it = sCache.find(distribution.signature())) != sCache.end()) {
        // Return cached result
        return it->second;

    } else {
        // Create a new cache
        AshikhminCache cache(distribution);
        std::pair<AshikhminCacheMap::iterator, bool> result =
                sCache.insert(AshikhminCacheMap::value_type(distribution.signature(), cache));
        Assert(result.second == true);
        fprintf(stderr, "Add %s to cache\n", result.first->first.c_str());  // TODO: test
        return result.first->second;
    }
}


Ashikhmin::Ashikhmin(const Spectrum &reflectance, Fresnel *f,
                       MicrofacetDistribution *d)   // TODO: use reference
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
     R(reflectance), mDistribution(d), fresnel(f),
     mCache(AshikhminCache::get(*d))
{
}

Spectrum
Ashikhmin::f(const Vector &woInput, const Vector &wiInput) const
{
    // In PBRT's implementation, woInput and wiInput is not guaranteed to be at the same side of the shading normal
    // My implementation of BlinnForAshikhmin::D() has domain for the whole sphere, and it's centered at
    // the local normal (0, 0, 1). Therefore, we need to inverse them if they are on the other side.
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

    Vector wh = wi + wo;    // Now wh is guaranteed to be at the same side as the local shading normal (0, 0, 1)
    if (wh.x == 0. && wh.y == 0. && wh.z == 0.) {
        return Spectrum(0.f);
    }
    wh = Normalize(wh);
    float i_dot_h = Dot(wi, wh);
    Spectrum F = fresnel->Evaluate(i_dot_h);
    //float avgNH = averageNH();  // TODO: no need to recompute everytime
    //float g_wi = gFactor(wi), g_wo = gFactor(wo);
    //fprintf(stderr, "%f %f %f\n", avgNH, g_wi, g_wo);
    float avgNH = 1;
    float g_wi = 1, g_wo = 1;
    // TODO: we need to make sure distribution->D(wh) actually returns a valid pdf; that is, it integrates to one over whole sphere
    return R * mDistribution->D(wh) * avgNH * F /
               (4.f * g_wi * g_wo);
}

Spectrum
Ashikhmin::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const 
{
    mDistribution->Sample_f(wo, wi, u1, u2, pdf);
    if (!SameHemisphere(wo, *wi)) {
        return Spectrum(0.f);
    }
    return f(wo, *wi);
}

float
Ashikhmin::Pdf(const Vector &wo, const Vector &wi) const 
{
    if (!SameHemisphere(wo, wi)) {
        return 0.f;
    }
    return mDistribution->Pdf(wo, wi);
}

float
Ashikhmin::averageNH(void) const
{
    // TODO
    return 1.f;
}

float
Ashikhmin::computeAverageNH(const MicrofacetDistribution &distribution)
{
    int nregions, neval, fail;
    double integral[1], error[1], prob[1];
    Cuhre(2, 1, averageNHIntegrand, (void*)&distribution, 1,
            1e-3, 1e-12, 0,
            0, 50000, 0,
            NULL,
            &nregions, &neval, &fail, integral, error, prob);

    // TODO: test
    //printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    //            nregions, neval, fail);
    //printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //            integral[0], error[0], prob[0]);

    return integral[0];
}

int
Ashikhmin::averageNHIntegrand(const int *ndim, const double xx[],
            const int *ncomp, double ff[], void *userdata)
{
    const MicrofacetDistribution *distribution =
                reinterpret_cast<const MicrofacetDistribution*>(userdata);
    float phi = xx[0] * 2*M_PI,
          theta = xx[1] * M_PI;

    const float costheta = cosf(theta), 
                sintheta = sinf(theta);

    Vector H = SphericalDirection(sintheta, costheta, phi);
    float NdotH = costheta;

    ff[0] = NdotH * distribution->D(H) * sintheta;
    ff[0] *= M_PI * (2*M_PI);

    return 0;
}

float
Ashikhmin::gFactor(const Vector &v) const
{
    // TODO
    return 1.f;
}

float
Ashikhmin::computeGFactor(const Vector &v, const MicrofacetDistribution &distribution)
{
    Quaternion q(Vector(0, 0, 1), v);
    GFactorIntegrandData data;

    data.distribution = &distribution;
    data.nToV = q.ToTransform();

    int nregions, neval, fail;
    double integral[1], error[1], prob[1];
    Cuhre(2, 1, gFactorIntegrand, (void*)&data, 1,
            1e-3, 1e-12, 0,
            0, 50000, 0,
            NULL,
            &nregions, &neval, &fail, integral, error, prob);

    return integral[0];
}

int
Ashikhmin::gFactorIntegrand(const int *ndim, const double xx[],
            const int *ncomp, double ff[], void *userdata)
{
    GFactorIntegrandData *data =
            reinterpret_cast<GFactorIntegrandData*>(userdata);
    float phi = xx[0] * 2*M_PI,
          theta = xx[1] * M_PI;

    // Find H in the space defined by V
    float costheta = cosf(theta), 
          sintheta = sinf(theta);

    Vector H = SphericalDirection(sintheta, costheta, phi);
    float VdotH = costheta;

    // Find H in the space defined by N
    Vector H_N = data->nToV(H);

    ff[0] = VdotH * data->distribution->D(H_N) * sintheta;
    ff[0] *= M_PI * (2*M_PI);

    return 0;
}

void
Ashikhmin::testSphVectorTransform(void)
{
    const int thetaRes = 1000, phiRes = 1000;
    float maxErrTheta = 0.f, maxErrPhi = 0.f;

    // The transform at boundary (theta = 0 or pi, phi = 0 or 2*pi) is not guaranteed to return the same value,
    // so they are skipped. Test result on 07/15/2014 is "maxErrTheta = 4.053116e-06, maxErrPhi = 4.768372e-07"
    for (int ti = 1; ti < thetaRes; ++ti) {
        float theta = M_PI * ti / thetaRes,
              sintheta = sinf(theta), 
              costheta = cosf(theta);
        for (int pi = 1; pi < phiRes; ++pi) {
            float phi = 2 * M_PI * pi / phiRes;
            Vector v = SphericalDirection(sintheta, costheta, phi);

            float errTheta = fabsf(SphericalTheta(v) - theta),
                  errPhi = fabsf(SphericalPhi(v) - phi);
            if (errTheta > maxErrTheta) {
                maxErrTheta = errTheta;
            }
            if (errPhi > maxErrPhi) {
                maxErrPhi = errPhi;
            }
        }
    }

    printf("Ashikhmin::testSphVectorTransform():\n"
           "  maxErrTheta = %e, maxErrPhi = %e\n", maxErrTheta, maxErrPhi);
}

void
Ashikhmin::testAverageNHAndFactor_g(void)
{
    const float BlinnExponent = 100.f;

    BlinnForAshikhmin distribution(BlinnExponent);
    float avgNH = Ashikhmin::computeAverageNH(distribution);
    float gFactorForN = Ashikhmin::computeGFactor(Vector(0, 0, 1), distribution),
          gFactorFor45deg = Ashikhmin::computeGFactor(Vector(1, 0, 1)/sqrtf(2), distribution);

    printf("For Blinn exponent %.2f:\n", BlinnExponent);
    printf("  Average dot(N,H) is %f\n", avgNH);
    printf("  g factor for N is %f\n", gFactorForN);
    printf("  g factor for 45 degree is %f\n", gFactorFor45deg);
}
