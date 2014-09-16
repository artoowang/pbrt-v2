
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

#include <map>

#include "reflection.h"
#include "montecarlo.h" // For Distribution1D
// TODO: remove this - not necessay anymore
//#include "mipmap.h"

// Boost
#include <boost/thread/mutex.hpp>
// In OSX, I need to use boost::math::isnan because PBRT's code can't find isnan once I include mutex.hpp
// TODO: not sure if this works on other platform
#include <boost/math/special_functions/fpclassify.hpp>
using boost::math::isnan;

using std::map;


class InterpolatedGrid
{
public:
    InterpolatedGrid();

    bool init(float x1, float x2, int width, float y1, float y2, int height,
            const vector<float> &samples);
    float eval(float x, float y) const;

    int getXResolution(void) const { return mWidth; }
    int getYResolution(void) const { return mHeight; }
    float getX1(void) const { return mX1; }
    float getX2(void) const { return mX2; }
    float getY1(void) const { return mY1; }
    float getY2(void) const { return mY2; }

protected:
    float mX1, mX2, mY1, mY2;
    int mWidth, mHeight;
    vector<float> mSamples;
};


class AshikhminDistribution {
public:
    AshikhminDistribution() { }
    virtual ~AshikhminDistribution() { }
    virtual float D(const Vector &wh) const = 0;
    virtual void Sample_f(const Vector &wo, Vector *wi,
                          float u1, float u2, float *pdf) const = 0;
    virtual float Pdf(const Vector &wo, const Vector &wi) const = 0;
    virtual string signature(void) const { return ""; }

    void printUnitTestResults(void) const;
    void writePdfImage(const Vector &wo, int thetaRes, int phiRes, const string &filepath) const;
    void writeDImage(int thetaRes, int phiRes, const string &filepath) const;

protected:
    // Following methods is for unit testing
    float integrateD(void) const;
    float integratePdf(const Vector &wo) const;

    static int DIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
                unsigned /*fdim*/, double *fval);
    static int PdfIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
                    unsigned /*fdim*/, double *fval);

    struct PdfIntegrandData
    {
        const AshikhminDistribution *distribution;
        Vector wo;
    };
};

// D is normalized for sphere (instead of the projected hemisphere as in Blinn)
// I.e., Integrate[D(wh), {wh in sphere}] = 1
// Currently Sample_f() and Pdf() implementation is the same as Blinn
// I think Pdf() is against wi instead of wh, but still need to confirm that
class BlinnForAshikhmin : public AshikhminDistribution
{
public:
    BlinnForAshikhmin(float e);
    virtual float D(const Vector &wh) const;
    virtual void Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wo, const Vector &wi) const;
    virtual string signature(void) const;

private:
    float exponent;
};

class TabulatedDistribution : public AshikhminDistribution
{
public:
    virtual ~TabulatedDistribution();
    virtual float D(const Vector &wh) const;
    virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wi, const Vector &wo) const;
    virtual string signature(void) const;

    void saveToFile(const string &filePath) const;

    static const TabulatedDistribution& get(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes);
    static const TabulatedDistribution& get(const string &filename);

    // TODO: for test
    static void test(void);

private:
    int mThetaRes, mPhiRes;
    vector<Distribution1D*> mThetaDists;
    Distribution1D *mPhiDist;
    vector<float> mPdf; // This could be combined with mPhiDist and mThetaDists if they allow access to the data
    InterpolatedGrid mData;
    string mSignature;

    // Object is only created by static member function
    TabulatedDistribution();
    // Do not allow copy
    TabulatedDistribution(const TabulatedDistribution&);
    TabulatedDistribution& operator=(const TabulatedDistribution&);

    void init(const vector<float> &dGrid);
    void initFromDistribution(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes);
    void initFromFile(const string &filePath);

    // Utility functions
    int getThetaRes(void) const { Assert(mThetaRes > 0); return mThetaRes; }
    int getPhiRes(void) const { Assert(mThetaRes > 0); return mPhiRes; }
    // The following two functions return the theta/phi value of the specified
    // cell. (Return the center of the cell by default.)
    float indexToTheta(int i, float offset = 0.5f) const { return M_PI * (i + offset) / getThetaRes(); }
    float indexToPhi(int i, float offset = 0.5f) const { return 2.f * M_PI * (i + offset) / getPhiRes(); }
    // The following two functions return the index of the cell where the
    // theta/phi value resides
    int thetaToIndex(float theta) const { return std::min(std::max((int)(theta/M_PI * getThetaRes()), 0), getThetaRes()-1); }
    int phiToIndex(float phi) const { return std::min(std::max((int)(phi/(2.f*M_PI) * getPhiRes()), 0), getPhiRes()-1); }
    // We use theta as the 1st dimension
    int gridIndex(int thetaIndex, int phiIndex) const { return phiIndex * getThetaRes() + thetaIndex; }

    // Static
    typedef map<string, TabulatedDistribution*> TabulatedDistributionMap;
    static TabulatedDistributionMap sCache;
    static boost::mutex sMutex;

    static string buildSignature(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes);
    static string buildSignature(const string &filePath);
};


class AshikhminCache
{
public:
    AshikhminCache();   // This is necessary for STL map
    virtual ~AshikhminCache();

    float gFactor(const Vector &v) const;
    float averageNH(void) const;

    static const AshikhminCache& get(const AshikhminDistribution &distribution);

    // TODO: test
    static int sGridResolution;

protected:
    float mAvgNH;
    //MIPMap<float> *mGFactorGrid;
    InterpolatedGrid mGFactorGrid;

    // A filled cache is only created by static member function
    AshikhminCache(const AshikhminDistribution &distribution);

    void initGGrid(int thetaRes, int phiRes, const AshikhminDistribution &distribution);

    typedef map<string, AshikhminCache*> AshikhminCacheMap;
    static AshikhminCacheMap sCache;
    static boost::mutex sMutex;

private:
    // This class does not allow copy
    AshikhminCache(const AshikhminCache&);
    AshikhminCache& operator=(const AshikhminCache&);
};

class Ashikhmin : public BxDF
{
public:
    Ashikhmin(const Spectrum &reflectance, Fresnel *f,
            const AshikhminDistribution &d);
    Spectrum f(const Vector &wo, const Vector &wi) const;
    Spectrum Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const;
    float Pdf(const Vector &wo, const Vector &wi) const;

    // TODO: for test
    bool mUseUniformSampling;

    // Utility functions
    static float computeGFactor(const Vector &v, const AshikhminDistribution &distribution);
    static float computeAverageNH(const AshikhminDistribution &distribution);

    // For tests
    static void testSphVectorTransform(void);
    static void testAverageNHAndFactor_g(void);

private:
    float averageNH(void) const;
    float gFactor(const Vector &v) const;

    static int averageNHIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
            unsigned /*fdim*/, double *fval);
    static int gFactorIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
            unsigned /*fdim*/, double *fval);

    Spectrum R;
    const AshikhminDistribution &mDistribution;
    Fresnel *fresnel;
    // Note: this needs to be reference, because if it's an object, it's destructor might not be called
    //       upon destruction if this object (Ashikhmin) is allocated by BSDF_ALLOC() (Need confirm)
    //       To avoid confusion, I have set its copy constructor to private
    const AshikhminCache &mCache;

    struct GFactorIntegrandData
    {
        const AshikhminDistribution *distribution;
        Vector v;
    };
};

#endif // PBRT_CORE_ASHIKHMIN_H
