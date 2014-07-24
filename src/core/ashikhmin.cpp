
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
#include "cubature-1.0/cubature.h"
#include "cuba.h"

using std::stringstream;

const double sCubatureRelError = 1e-4;

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
boost::mutex AshikhminCache::sMutex;

AshikhminCache::AshikhminCache() :
        mAvgNH(-1.f), mGFactorGrid(NULL)
{
}

AshikhminCache::AshikhminCache(const MicrofacetDistribution &distribution) :
        mAvgNH(Ashikhmin::computeAverageNH(distribution)), mGFactorGrid(NULL)
{
    // TODO: test
    fprintf(stderr, "Cache for %s created.\n", distribution.signature().c_str());

    // Initialize grid of factor g
    // TODO: currently using 32x32
    const int thetaRes = 32, phiRes = 32;
    initGGrid(thetaRes, phiRes, distribution);
}

AshikhminCache::~AshikhminCache()
{
    if (mGFactorGrid != NULL) {
        delete mGFactorGrid;
        mGFactorGrid = NULL;
    }
}

void
AshikhminCache::initGGrid(int thetaRes, int phiRes, const MicrofacetDistribution &distribution)
{
    float *gGrid = new float[thetaRes * phiRes];

    Assert(gGrid != NULL);  // TODO: error handling

    // The center of the upper-left pixel (x,y) = (0,0) is mapped to
    // (s,t) = (0,0)+half_pixel_size, while the lower-right pixel (x,y) = (width-1,height-1)
    // is mapped to (s,t) = (1,1)-half_pixel_size
    // I.e., s = (x+0.5)/width, t = (y+0.5)/height
    // Then, theta = s*pi, phi = t*2*pi
    for (int x = 0; x < thetaRes; ++x) {
        const float s = (x + 0.5f) / thetaRes,
                    theta = M_PI * s,
                    costheta = cosf(theta),
                    sintheta = sinf(theta);
        for (int y = 0; y < phiRes; ++y) {
            const float t = (y + 0.5f) / phiRes,
                        phi = 2.f * M_PI * t;

            Vector v = SphericalDirection(sintheta, costheta, phi);
            gGrid[y*thetaRes + x] = Ashikhmin::computeGFactor(v, distribution);
        }
    }

    mGFactorGrid = new MIPMap<float>(thetaRes, phiRes, gGrid, false, 8.f, TEXTURE_CLAMP);
    delete [] gGrid;
}

float
AshikhminCache::gFactor(const Vector &v) const
{
    Assert(mGFactorGrid != NULL);
    const float theta = SphericalTheta(v),
                s = theta * INV_PI,
                phi = SphericalPhi(v),
                t = phi * INV_TWOPI;
    return mGFactorGrid->Lookup(s, t);
}

float
AshikhminCache::averageNH(void) const
{
    return mAvgNH;
}

const AshikhminCache&
AshikhminCache::get(const MicrofacetDistribution &distribution)
{
    boost::mutex::scoped_lock scopedLock(sMutex);
    AshikhminCacheMap::const_iterator it;
    if ((it = sCache.find(distribution.signature())) != sCache.end()) {
        // Return cached result
        Assert(it->second != NULL);
        return *(it->second);

    } else {
        // Create a new cache
        AshikhminCache *cache = new AshikhminCache(distribution);   // TODO: is there a way to release them upon program exit?
        Assert(cache != NULL); // TODO: error handling
        std::pair<AshikhminCacheMap::iterator, bool> result =
                sCache.insert(AshikhminCacheMap::value_type(distribution.signature(), cache));
        Assert(result.second == true);
        Assert(result.first->second != NULL);
        return *(result.first->second);
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
    float avgNH = averageNH();
    float g_wi = gFactor(wi), g_wo = gFactor(wo);
    // TODO: remove
    //fprintf(stderr, "%f %f %f\n", avgNH, g_wi, g_wo);
    //float avgNH = 1;
    //float g_wi = 1, g_wo = 1;
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
    return mCache.averageNH();
}

float
Ashikhmin::computeAverageNH(const MicrofacetDistribution &distribution)
{
    // First dimension is phi, second theta
    const double xmin[2] = {0, 0},
                 xmax[2] = {2*M_PI, M_PI};
    double val = 0.f, error = 0.f;
    int ret = hcubature(1, averageNHIntegrand, (void*)&distribution,
            2, xmin, xmax,
            0, 0, sCubatureRelError,
            ERROR_INDIVIDUAL, &val, &error);

    Assert(ret == 0);

    // TODO: test
    //printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    //            nregions, neval, fail);
    //printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //            integral[0], error[0], prob[0]);

    return val;
}

// ndim should be 2, fdim should be 1
int
Ashikhmin::averageNHIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
        unsigned /*fdim*/, double *fval)
{
    const MicrofacetDistribution *distribution =
                reinterpret_cast<const MicrofacetDistribution*>(fdata);
    const float phi = x[0], theta = x[1];
    const float costheta = cosf(theta), 
                sintheta = sinf(theta);

    Vector H = SphericalDirection(sintheta, costheta, phi);
    float NdotH = costheta;

    fval[0] = NdotH * distribution->D(H) * sintheta;

    return 0;
}

float
Ashikhmin::gFactor(const Vector &v) const
{
    return mCache.gFactor(v);
}

float
Ashikhmin::computeGFactor(const Vector &v, const MicrofacetDistribution &distribution)
{
    Quaternion q(Vector(0, 0, 1), v);
    GFactorIntegrandData data;

    data.distribution = &distribution;
    data.nToV = q.ToTransform();

    // First dimension is phi, second theta
    const double xmin[2] = {0, 0},
                 xmax[2] = {2*M_PI, M_PI};
    double val = 0.f, error = 0.f;
    int ret = hcubature(1, gFactorIntegrand, (void*)&data,
            2, xmin, xmax,
            0, 0, sCubatureRelError,
            ERROR_INDIVIDUAL, &val, &error);

    Assert(ret == 0);

    // TODO: test
    //printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    //            nregions, neval, fail);
    //printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
    //            integral[0], error[0], prob[0]);

    return val;
}

// ndim should be 2, fdim should be 1
int
Ashikhmin::gFactorIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
        unsigned /*fdim*/, double *fval)
{
    GFactorIntegrandData *data =
            reinterpret_cast<GFactorIntegrandData*>(fdata);
    const float phi = x[0], theta = x[1];

    // Find H in the space defined by V
    const float costheta = cosf(theta),
          sintheta = sinf(theta);

    Vector H = SphericalDirection(sintheta, costheta, phi);
    float VdotH = costheta;

    // Find H in the space defined by N
    Vector H_N = data->nToV(H);

    fval[0] = VdotH * data->distribution->D(H_N) * sintheta;

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

    fprintf(stderr, "Ashikhmin::testSphVectorTransform():\n"
                    "  maxErrTheta = %e, maxErrPhi = %e\n", maxErrTheta, maxErrPhi);
}

void
Ashikhmin::testAverageNHAndFactor_g(void)
{
    const float BlinnExponent = 100.f;
    const Vector N(0, 0, 1), V45(1/sqrtf(2), 0, 1/sqrtf(2));

    BlinnForAshikhmin distribution(BlinnExponent);
    float avgNH = Ashikhmin::computeAverageNH(distribution);
    float gFactorForN = Ashikhmin::computeGFactor(N, distribution),
          gFactorFor45deg = Ashikhmin::computeGFactor(V45, distribution);

    // Test cache
    const AshikhminCache &cache = AshikhminCache::get(distribution);
    float cacheAvgNH = cache.averageNH();
    float cacheGFactorForN = cache.gFactor(N),
          cacheGFactorFor45deg = cache.gFactor(V45);

    fprintf(stderr, "For Blinn exponent %.2f:\n", BlinnExponent);
    fprintf(stderr, "  Average dot(N,H) is %f (%f)\n", avgNH, cacheAvgNH);
    fprintf(stderr, "  g factor for N is %f (%f)\n", gFactorForN, cacheGFactorForN);
    fprintf(stderr, "  g factor for 45 degree is %f (%f)\n", gFactorFor45deg, cacheGFactorFor45deg);
}
