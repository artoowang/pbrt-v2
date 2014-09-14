
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

// TODO: remove if not used
#include <tiffio.h>     // For TIFF output

// For EXR support
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>

// For numerical integration
#include "cubature-1.0/cubature.h"

#include <boost/format.hpp>

using std::stringstream;

// For EXR support
using namespace Imf;
using namespace Imath;

const double sCubatureRelError = 1e-4;
const double sCubatureAbsError = 1e-9;  // This is necessary for integral that leads to zero
const float sSmallValue = 1e-6f;
const bool sPrintGGrid = false;

// ----------------------------------------------------------------------------

// Copied from tifftoexr.cpp
void WriteEXR(const char *name, float *frgba, int xRes, int yRes, bool hasAlpha)
{
    Header header(xRes, yRes);
    header.channels().insert("R", Channel (HALF));
    header.channels().insert("G", Channel (HALF));
    header.channels().insert("B", Channel (HALF));
    if (hasAlpha)
    header.channels().insert("A", Channel (HALF));
    int stride = hasAlpha ? 4 : 3;

    half *rgba = new half[xRes*yRes * stride];
    for (int i = 0; i < xRes*yRes * stride; ++i)
    rgba[i] = frgba[i];

    FrameBuffer fb;
    fb.insert("R", Slice(HALF, (char *)rgba, stride*sizeof(half),
             stride*xRes*sizeof(half)));
    fb.insert("G", Slice(HALF, (char *)rgba+sizeof(half), stride*sizeof(half),
             stride*xRes*sizeof(half)));
    fb.insert("B", Slice(HALF, (char *)rgba+2*sizeof(half), stride*sizeof(half),
             stride*xRes*sizeof(half)));
    if (hasAlpha)
    fb.insert("A", Slice(HALF, (char *)rgba+3*sizeof(half), stride*sizeof(half),
                 stride*xRes*sizeof(half)));

    OutputFile file(name, header);
    file.setFrameBuffer(fb);
    file.writePixels(yRes);
}

// rgba: [0, 255]
void WriteTIFF(const char *name, float *rgba, int XRes, int YRes, bool hasAlpha)
{
    // Open 8-bit TIFF file for writing
    TIFF *tiff = TIFFOpen(name, "w");
    if (!tiff) {
    fprintf(stderr, "Unable to open TIFF %s for writing", name);
    return;
    }

    int nChannels = hasAlpha ? 4 : 3;
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, nChannels);
    if (hasAlpha) {
    short int extra[] = { EXTRASAMPLE_ASSOCALPHA };
    TIFFSetField(tiff, TIFFTAG_EXTRASAMPLES, (short)1, extra);
    }
    // Write image resolution information
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, XRes);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, YRes);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    // Set Generic TIFF Fields
    TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(tiff, TIFFTAG_XRESOLUTION, 1.f);
    TIFFSetField(tiff, TIFFTAG_YRESOLUTION, 1.f);
    TIFFSetField(tiff, TIFFTAG_RESOLUTIONUNIT, (short)1);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
    // Write 8-bit scanlines
    unsigned char *buf = new unsigned char[nChannels * XRes];
    for (int y = 0; y < YRes; ++y) {
    unsigned char *bufp = buf;
    for (int x = 0; x < XRes; ++x) {
        // Pack 8-bit pixels samples into buf
        for (int s = 0; s < nChannels; ++s)
        *bufp++ = (unsigned char)*rgba++;
    }
    TIFFWriteScanline(tiff, buf, y, 1);
    }
    // Close 8-bit TIFF file
    delete[] buf;
    TIFFClose(tiff);
}

// ----------------------------------------------------------------------------

InterpolatedGrid::InterpolatedGrid() :
        mX1(0), mX2(0), mY1(0), mY2(0),
        mWidth(-1), mHeight(-1)
{
}

// samples is assumed to be formatted in row-major order
// i.e., samples[x][y] = samples[y*width+x]
bool
InterpolatedGrid::init(float x1, float x2, int width,
        float y1, float y2, int height,
        const vector<float> &samples)
{
    if (x2 < x1 || y2 < y1 || width <= 0 || height <= 0) {
        return false;
    }

    if (samples.size() < static_cast<vector<float>::size_type>(width * height)) {
        return false;
    }

    mX1 = x1; mX2 = x2; mY1 = y1; mY2 = y2;
    mWidth = width; mHeight = height;
    mSamples.assign(samples.begin(), samples.begin() + mWidth*mHeight);

    return true;
}

// (x, y) is in the coordinates defined by [mX1, mX2] x [mY1, mY2]
float
InterpolatedGrid::eval(float x, float y) const
{
    Assert(mWidth > 0 && mHeight > 0);

    // Transform x, y into pixel coordinate system:
    //   the upper-left samples is at (0, 0), lower-right is at
    //   (mWidth-1, mHeight-1), and the upper-left corner is
    //   (-0.5, -0.5), the lower-right corner (mWidth-0.5, mHeight-0.5)
    x = (x - mX1) / (mX2 - mX1) * mWidth - 0.5f;
    y = (y - mY1) / (mY2 - mY1) * mHeight - 0.5f;
    int x1 = (int)x,
        y1 = (int)y,
        x2 = x1 + 1,
        y2 = y1 + 1;
    float dx = x - x1,
          dy = y - y1;
    x1 = std::max(0, x1);
    y1 = std::max(0, y1);
    x2 = std::min(mWidth-1, x2);
    y2 = std::min(mHeight-1, y2);
    return (1-dx) * (1-dy) * mSamples[mWidth*y1 + x1]
         + dx     * (1-dy) * mSamples[mWidth*y1 + x2]
         + (1-dx) * dy     * mSamples[mWidth*y2 + x1]
         + dx     * dy     * mSamples[mWidth*y2 + x2];
}

// ----------------------------------------------------------------------------

void
AshikhminDistribution::printUnitTestResults(void) const
{
    fprintf(stderr, "Running unit tests for %s\n", signature().c_str());
    fprintf(stderr, "  Integrate D(wh) over wh (should be one): %f\n",
            integrateD());
    // Currently just use some arbitrary wo. Note the (0, 0, 1) is the
    // local surface normal
    const Vector wo = Normalize(Vector(1, -2, 3));
    fprintf(stderr, "  Using wo = (%.2f, %.2f, %.2f)\n", wo.x, wo.y, wo.z);
    fprintf(stderr, "  Integrate Pdf(wo, wi) over wi (should be one): %f\n",
            integratePdf(wo));

    // Output PDF/D images
    writePdfImage(wo, 512, 512, signature() + ".pdf.exr");
    writeDImage(512, 512, signature() + ".D.exr");
}

void
AshikhminDistribution::writePdfImage(const Vector &wo, int thetaRes, int phiRes, const string &filepath) const
{
    vector<float> data(thetaRes * phiRes, 0.f);

    // Fill in the PDF values
    for (int y = 0; y < thetaRes; ++y) {
        const float theta = M_PI * (y + 0.5f) / thetaRes,
                                costheta = cosf(theta),
                                sintheta = sinf(theta);

        for (int x = 0; x < phiRes; ++x) {
            const float phi = 2.f * M_PI * (x + 0.5f) / phiRes;
            const Vector &wi = SphericalDirection(sintheta, costheta, phi);
            data[y*thetaRes + x] = Pdf(wo, wi);
        }
    }

    // Convert to 3 channels
    vector<float> pixels(thetaRes * phiRes * 3, 0.f);
    for (int i = 0; i < thetaRes * phiRes; ++i) {
        pixels[i*3] = pixels[i*3+1] = pixels[i*3+2] = data[i];
    }

    WriteEXR(filepath.c_str(), pixels.data(), phiRes, thetaRes, false);
}

void
AshikhminDistribution::writeDImage(int thetaRes, int phiRes, const string &filepath) const
{
    vector<float> data(thetaRes * phiRes, 0.f);

    // Fill in the D values
    for (int y = 0; y < thetaRes; ++y) {
        const float theta = M_PI * (y + 0.5f) / thetaRes,
                                costheta = cosf(theta),
                                sintheta = sinf(theta);

        for (int x = 0; x < phiRes; ++x) {
            const float phi = 2.f * M_PI * (x + 0.5f) / phiRes;
            const Vector &wh = SphericalDirection(sintheta, costheta, phi);
            data[y*thetaRes + x] = D(wh);
        }
    }

    // Convert to 3 channels
    vector<float> pixels(thetaRes * phiRes * 3, 0.f);
    for (int i = 0; i < thetaRes * phiRes; ++i) {
        pixels[i*3] = pixels[i*3+1] = pixels[i*3+2] = data[i];
    }

    WriteEXR(filepath.c_str(), pixels.data(), phiRes, thetaRes, false);
}

float
AshikhminDistribution::integrateD(void) const
{
    // First dimension is phi, second theta
    const double xmin[2] = {0, 0},
                   xmax[2] = {2.f*M_PI, M_PI};
    double val = 0.f, error = 0.f;
#ifdef DEBUG
    int ret =
#endif
    hcubature(1, DIntegrand, (void*)this,
            2, xmin, xmax,
            0, 0, sCubatureRelError,
            ERROR_INDIVIDUAL, &val, &error);

    Assert(ret == 0);

    return val;
}

// ndim should be 2, fdim should be 1
int
AshikhminDistribution::DIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
        unsigned /*fdim*/, double *fval)
{
    const AshikhminDistribution *distribution =
                reinterpret_cast<const AshikhminDistribution*>(fdata);
    const float phi = x[0], theta = x[1];
    const float costheta = cosf(theta),
                 sintheta = sinf(theta);

    Vector H = SphericalDirection(sintheta, costheta, phi);
    fval[0] = distribution->D(H) * sintheta;

    return 0;
}

float
AshikhminDistribution::integratePdf(const Vector &wo) const
{
    PdfIntegrandData data;
    data.distribution = this;
    data.wo = wo;

    // First dimension is phi, second theta
    const double xmin[2] = {0, 0},
                   xmax[2] = {2.f*M_PI, M_PI};
    double val = 0.f, error = 0.f;
#ifdef DEBUG
    int ret =
#endif
    hcubature(1, PdfIntegrand, (void*)&data,
            2, xmin, xmax,
            0, 0, sCubatureRelError,
            ERROR_INDIVIDUAL, &val, &error);

    Assert(ret == 0);

    return val;
}

// ndim should be 2, fdim should be 1
int
AshikhminDistribution::PdfIntegrand(unsigned /*ndim*/, const double *x, void *fdata,
        unsigned /*fdim*/, double *fval)
{
    const PdfIntegrandData *data =
                reinterpret_cast<const PdfIntegrandData*>(fdata);
    const float phi = x[0], theta = x[1];
    const float costheta = cosf(theta),
                 sintheta = sinf(theta);

    Vector wi = SphericalDirection(sintheta, costheta, phi);
    fval[0] = data->distribution->Pdf(data->wo, wi) * sintheta;

    return 0;
}

// ----------------------------------------------------------------------------

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
    // Note: unlike class Blinn, we want to make sure there is no microfacet facing downward.
    //       This makes Integrate[D(wh), {wh in sphere}] = 1
    float costhetah = max(CosTheta(wh), 0.f);
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
    // Note: not sure why we don't check this earlier. We need to do this because
    //       otherwise the PDF becomes negative. When dot(wo, wh) <= 0, it means
    //       the angle between wo and wh is larger than 90 degree, and therefore
    //       wi and wo won't be in the same hemisphere anyway. (For normal BRDFs,
    //       this usually means the contribution is zero.)
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

// ----------------------------------------------------------------------------
// TODO: remove

TabulatedDistributionTest::TabulatedDistributionMap TabulatedDistributionTest::sCache;
boost::mutex TabulatedDistributionTest::sMutex;

TabulatedDistributionTest::TabulatedDistributionTest() :
        mDistribution(NULL), mSignature("TabulatedDistributionTest::uninitialized")
{
}

// TODO: currently, if this is allocated by BSDF_ALLOC or TabulatedDistributionTest::get(),
//       it won't be destructed. (For the former, it's because the memory blocks are
//       released directly by MemoryAreana; for the later, it's because the static cache
//       currently doesn't release itself at the end of the application
TabulatedDistributionTest::~TabulatedDistributionTest()
{
    if (mDistribution != NULL) {
        delete mDistribution;
    }
}

void
TabulatedDistributionTest::initFromDistribution(const AshikhminDistribution& srcDistribution, int thetaRes, int phiRes)
{
    // We parse through the source distribution by a thetaRes x phiRes grid,
    // over [0,pi] x [0,2*pi]. The samples are sampled at the center of each cell
    // - data: the PDF (against half vector wh over sphere) returned by srcDistribution.D()
    // - pdf: the PDF against 2D pair (theta, phi) over [0,pi] x [0,2*pi], which is
    //        D(h) * sin(theta). This is the value used to drive Distribution1D
    vector<float> data(thetaRes * phiRes), pdf(thetaRes * phiRes);
    for (int x = 0; x < thetaRes; ++x) {
        const float s = (x + 0.5f) / thetaRes,
                    theta = M_PI * s,
                    costheta = cosf(theta),
                    sintheta = sinf(theta);

        for (int y = 0; y < phiRes; ++y) {
            const float t = (y + 0.5f) / phiRes,
                        phi = 2.f * M_PI * t;
            const Vector &h = SphericalDirection(sintheta, costheta, phi);
            const float D = srcDistribution.D(h);
            data[y*thetaRes + x] = D;
            pdf[y*thetaRes + x] = D * sintheta;
        }
    }

    // The 1D distribution to select cells against a PDF in (theta, phi) domain
    // TODO: one problem of using 1D domain is that we might lose the stratified
    //       property of our random numbers
    if (mDistribution != NULL) {
        delete mDistribution;
    }
    mDistribution = new Distribution1D(pdf.data(), thetaRes * phiRes);
    Assert(mDistribution != NULL);

    // The actual PDF of wh in sphere domain
    mData.init(0, M_PI, thetaRes, 0, 2*M_PI, phiRes, data);

    mSignature = buildSignature(srcDistribution, thetaRes, phiRes);

    // TODO: test
    fprintf(stderr, "%s created.\n", mSignature.c_str());
}

float
TabulatedDistributionTest::D(const Vector &wh) const
{
    const float theta = SphericalTheta(wh),
                phi = SphericalPhi(wh);
    return mData.eval(theta, phi);
}

void
TabulatedDistributionTest::Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const
{
    // It seems we always call Sample_f() with a pdf, and we use it to
    // determine if the sample is valid
    Assert(pdf != NULL);

    const int thetaRes = mData.getXResolution(),
              phiRes = mData.getYResolution(),
              index = mDistribution->SampleDiscrete(u1, NULL),
              thetaIndex = index % thetaRes,
              phiIndex = index / thetaRes;

    // TODO: we should jitter the wi within the selected cell
    const float theta = ((thetaIndex + 0.5f) / thetaRes) * M_PI,
                phi = ((phiIndex + 0.5f) / phiRes) * 2.f * M_PI,
                sinTheta = sinf(theta), cosTheta = cosf(theta);
    if (theta > 0.5f * M_PI) {
        fprintf(stderr, "theta = %f\n", theta);
    }
    const Vector wh = SphericalDirection(sinTheta, cosTheta, phi);
    const float dotHO = Dot(wh, wo);
    if (dotHO < 0) {
        //fprintf(stderr, "dotHO = %f\n", dotHO);
        *pdf = 0.f;
        return;
    }

    // Compute incident direction by reflecting about wh
    *wi = -wo + 2.f * dotHO * wh;
    *pdf = Pdf(wo, *wi);    // TODO: can be optimized
}

float
TabulatedDistributionTest::Pdf(const Vector &wo, const Vector &wi) const
{
    const Vector wh = Normalize(wo + wi);
    const float dotHO = Dot(wh, wo);
    if (dotHO > 0.f) {
        return D(wh) / (4.f * dotHO);
    } else {
        // It should be very unlikely we are here, but let's check anyway
        return 0.f;
    }
}

string
TabulatedDistributionTest::signature(void) const
{
    return mSignature;
}

string
TabulatedDistributionTest::buildSignature(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes)
{
    stringstream ss;
    ss << "TabulatedDistributionTest:"
            << srcDistribution.signature() << ":"
            << thetaRes << "x" << phiRes;
    return ss.str();
}

const TabulatedDistributionTest&
TabulatedDistributionTest::get(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes)
{
    boost::mutex::scoped_lock scopedLock(sMutex);
    TabulatedDistributionMap::const_iterator it;
    const string &signature = buildSignature(srcDistribution, thetaRes, phiRes);

    if ((it = sCache.find(signature)) != sCache.end()) {
        // Return cached result
        Assert(it->second != NULL);
        return *(it->second);

    } else {
        // Create a new TabulatedDistributionTest
        TabulatedDistributionTest *distribution = new TabulatedDistributionTest;   // TODO: is there a way to release them upon program exit?
        Assert(distribution != NULL); // TODO: error handling
        distribution->initFromDistribution(srcDistribution, thetaRes, phiRes);

        std::pair<TabulatedDistributionMap::iterator, bool> result =
                sCache.insert(TabulatedDistributionMap::value_type(distribution->signature(), distribution));
        Assert(result.second == true);
        Assert(result.first->second != NULL);
        Assert(distribution->signature() == signature);

        return *(result.first->second);
    }
}

// ----------------------------------------------------------------------------

TabulatedDistribution::TabulatedDistributionMap TabulatedDistribution::sCache;
boost::mutex TabulatedDistribution::sMutex;

TabulatedDistribution::TabulatedDistribution() :
        mPhiDist(NULL), mSignature("TabulatedDistribution::uninitialized")
{
}

// TODO: currently, if this is allocated by BSDF_ALLOC or TabulatedDistribution::get(),
//       it won't be destructed. (For the former, it's because the memory blocks are
//       released directly by MemoryAreana; for the later, it's because the static cache
//       currently doesn't release itself at the end of the application
TabulatedDistribution::~TabulatedDistribution()
{
    if (mPhiDist != NULL) {
        delete mPhiDist;
    }
    for (vector<Distribution1D*>::iterator it = mThetaDists.begin();
            it != mThetaDists.end();
            ++it) {
        delete *it;
    }
}

void
TabulatedDistribution::initFromDistribution(const AshikhminDistribution& srcDistribution, int thetaRes, int phiRes)
{
    // We parse through the source distribution by a thetaRes x phiRes grid,
    // over [0,pi] x [0,2*pi]. The samples are sampled at the center of each cell
    // - data: the PDF (against half vector wh over sphere) returned by srcDistribution.D()
    // - mPdf: the PDF against 2D pair (theta, phi) over [0,pi] x [0,2*pi], which is
    //         D(h) * sin(theta). This is the value used to drive Distribution1D
    vector<float> data(thetaRes * phiRes), pdfy(phiRes);
    mPdf.resize(thetaRes * phiRes, 0.f);
    for (int y = 0; y < phiRes; ++y) {
        const float phi = 2.f * M_PI * (y + 0.5f) / phiRes;
        pdfy[y] = 0.f;

        for (int x = 0; x < thetaRes; ++x) {
            const float theta = M_PI * (x + 0.5f) / thetaRes,
                        costheta = cosf(theta),
                        sintheta = sinf(theta);
            const Vector &h = SphericalDirection(sintheta, costheta, phi);
            const float D = srcDistribution.D(h);

            data[y*thetaRes + x] = D;
            float pdf = D * sintheta;
            mPdf[y*thetaRes + x] = pdf;
            pdfy[y] += pdf;
        }
    }

    // The 1D distribution to select cells against a PDF in (theta, phi) domain
    // We first sample against phi using pdfy, then sample against theta using pdf+y*thetaRes
    // This utilizes both 2D samples (u1 and u2) and works with Halton sequence and
    // lowdiscrepency sampler. See Sample_f()
    Assert(mThetaDists.size() == 0 && mPhiDist == NULL);

    mPhiDist = new Distribution1D(pdfy.data(), phiRes);
    Assert(mPhiDist != NULL);

    mThetaDists.resize(phiRes, NULL);
    for (int y = 0; y < phiRes; ++y) {
        Distribution1D *dist = new Distribution1D(mPdf.data() + y*thetaRes, thetaRes);
        Assert(dist != NULL);
        mThetaDists[y] = dist;
    }

    // The actual PDF of wh in sphere domain
    mData.init(0, M_PI, thetaRes, 0, 2*M_PI, phiRes, data);

    mSignature = buildSignature(srcDistribution, thetaRes, phiRes);

    // TODO: test
    fprintf(stderr, "%s created.\n", mSignature.c_str());
    float pdfInt = 0.f;
    for (int i = 0; i < thetaRes * phiRes; ++i) {
        pdfInt += mPdf[i];
    }
    pdfInt *= M_PI / thetaRes * 2.f * M_PI / phiRes;
    fprintf(stderr, "PDF integrates to %f\n", pdfInt);
}

float
TabulatedDistribution::D(const Vector &wh) const
{
    const float theta = SphericalTheta(wh),
                phi = SphericalPhi(wh);
    return mData.eval(theta, phi);
}

// Notice the pdf of Sample_f() and Pdf() is w.r.t wi, not wh
// TODO: unlike existing PBRT implementations, we assume the normal of the
//       tangent frame is facing toward wo. (See Ashikhmin::f()). Therefore
//       we DON'T reverse it when wh is at the opposite side of wo, like
//       Blinn::Sample_f(); instead, we just check if Dot(wh, wo) < 0 - if it
//       is negative, the viewer won't see the facet and the sample is invalid
void
TabulatedDistribution::Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const
{
    /*// TODO: test "analytical" Sample_f() for BlinnForAshikhmin with exponent = 12.5
    // Compute sampled half-angle vector $\wh$ for Blinn distribution
    float exponent = 12.5f;
    float costheta = powf(u1, 1.f / (exponent+1));
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    float phi = u2 * 2.f * M_PI;
    Vector wh = SphericalDirection(sintheta, costheta, phi);
    if (!SameHemisphere(wo, wh)) {
        fprintf(stderr, "Reverse wh\n");
        wh = -wh;
    }

    if (Dot(wo, wh) < 0) {
        fprintf(stderr, "HO = %f\n", Dot(wo, wh));
        *pdf = 0.f;
        return;
    }

    // Compute incident direction by reflecting about $\wh$
    *wi = -wo + 2.f * Dot(wo, wh) * wh;
    *pdf = Pdf(wo, *wi);*/

    // It seems we always call Sample_f() with a pdf, and we use it to
    // determine if the sample is valid
    const int thetaRes = mData.getXResolution(),
              phiRes = mData.getYResolution();

    Assert(pdf != NULL && mPhiDist != NULL && mThetaDists.size() == phiRes);

    float phiU = 0.5f, thetaU = 0.5f;   // Used to jitter the sample inside the cells
    const int phiIndex = mPhiDist->SampleDiscrete(u1, NULL, &phiU),
              thetaIndex = mThetaDists[phiIndex]->SampleDiscrete(u2, NULL, &thetaU);

    const float theta = ((thetaIndex + thetaU) / thetaRes) * M_PI,
                phi = ((phiIndex + phiU) / phiRes) * 2.f * M_PI,
                sinTheta = sinf(theta), cosTheta = cosf(theta);
    const Vector wh = SphericalDirection(sinTheta, cosTheta, phi);
    const float dotHO = Dot(wh, wo);
    if (dotHO < 0) {
        *pdf = 0.f;
        return;
    }

    // Compute incident direction by reflecting about wh
    *wi = -wo + 2.f * dotHO * wh;
    *pdf = Pdf(wo, *wi);    // TODO: can be optimized
}

float
TabulatedDistribution::Pdf(const Vector &wo, const Vector &wi) const
{
    const Vector wh = Normalize(wo + wi);
    const float dotHO = Dot(wh, wo);

    if (dotHO > 0.f) {
        const float theta = SphericalTheta(wh),
                    phi = SphericalPhi(wh);
        const int thetaRes = mData.getXResolution(),
                  phiRes = mData.getYResolution();
        const int x = std::min(std::max((int)(theta/M_PI * thetaRes), 0), thetaRes-1),
                  y = std::min(std::max((int)(phi/(2.f*M_PI) * phiRes), 0), phiRes-1);
        // Notice we don't use the sin(theta) to transform PDF. Instead we use
        // the theta at the center of the cell. This also avoids sin(theta) = 0
        const float sintheta = sinf(M_PI * (x + 0.5f) / thetaRes);

        // Transform pdf against (theta, phi) space to solid angle space
        const float Dval = mPdf[y*thetaRes + x] / sintheta;
        // Transform pdf against wh to wi
        return Dval / (4.f * dotHO);

        // TODO: old computations
        //fprintf(stderr, "%f %f\n", D(wh), Dval);
        //return D(wh) / (4.f * dotHO);

    } else {
        // It should be very unlikely we are here, but let's check anyway
        return 0.f;
    }
}

string
TabulatedDistribution::signature(void) const
{
    return mSignature;
}

string
TabulatedDistribution::buildSignature(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes)
{
    stringstream ss;
    ss << "TabulatedDistribution:"
            << srcDistribution.signature() << ":"
            << thetaRes << "x" << phiRes;
    return ss.str();
}

const TabulatedDistribution&
TabulatedDistribution::get(const AshikhminDistribution &srcDistribution, int thetaRes, int phiRes)
{
    boost::mutex::scoped_lock scopedLock(sMutex);
    TabulatedDistributionMap::const_iterator it;
    const string &signature = buildSignature(srcDistribution, thetaRes, phiRes);

    if ((it = sCache.find(signature)) != sCache.end()) {
        // Return cached result
        Assert(it->second != NULL);
        return *(it->second);

    } else {
        // Create a new TabulatedDistribution
        TabulatedDistribution *distribution = new TabulatedDistribution;   // TODO: is there a way to release them upon program exit?
        Assert(distribution != NULL); // TODO: error handling
        distribution->initFromDistribution(srcDistribution, thetaRes, phiRes);

        // TODO: run unit tests
        srcDistribution.printUnitTestResults();
        distribution->printUnitTestResults();

        std::pair<TabulatedDistributionMap::iterator, bool> result =
                sCache.insert(TabulatedDistributionMap::value_type(distribution->signature(), distribution));
        Assert(result.second == true);
        Assert(result.first->second != NULL);
        Assert(distribution->signature() == signature);

        return *(result.first->second);
    }
}

void
TabulatedDistribution::test(void)
{
    fprintf(stderr, "Test Distribution1D:\n");
    float f[] = {1, 2, 3, 4};
    Distribution1D distribution(f, 4);
    fprintf(stderr, "  SampleDiscrete:\n");
    for (int i = 0; i < 10; ++i) {
        float pdf = 0.f, u = i / 10.f;
        int s = distribution.SampleDiscrete(u, &pdf);
        fprintf(stderr, "    u = %.2f, select %d, pdf = %.2f\n", u, s, pdf);
    }
    fprintf(stderr, "  SampleContinuous:\n");
    for (int i = 0; i < 10; ++i) {
        float pdf = 0.f, u = i / 10.f;
        float s = distribution.SampleContinuous(u, &pdf);
        fprintf(stderr, "    u = %.2f, select %.2f, pdf = %.2f\n", u, s, pdf);
    }
}

// ----------------------------------------------------------------------------

AshikhminCache::AshikhminCacheMap AshikhminCache::sCache;
boost::mutex AshikhminCache::sMutex;
int AshikhminCache::sGridResolution = 32;

AshikhminCache::AshikhminCache() :
        mAvgNH(-1.f)//, mGFactorGrid(NULL)
{
}

AshikhminCache::AshikhminCache(const AshikhminDistribution &distribution) :
        mAvgNH(Ashikhmin::computeAverageNH(distribution))//, mGFactorGrid(NULL)
{
    // TODO: test
    fprintf(stderr, "Cache for %s created.\n", distribution.signature().c_str());

    // Initialize grid of factor g
    int thetaRes, phiRes;
    thetaRes = phiRes = sGridResolution;   // TODO: test
    initGGrid(thetaRes, phiRes, distribution);
}

AshikhminCache::~AshikhminCache()
{
    //if (mGFactorGrid != NULL) {
    //    delete mGFactorGrid;
    //    mGFactorGrid = NULL;
    //}
}

void
AshikhminCache::initGGrid(int thetaRes, int phiRes, const AshikhminDistribution &distribution)
{
    vector<float> gGrid(thetaRes * phiRes);

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
        // TODO: test
        //fprintf(stderr, "theta = %f (%d/%d)\n", theta, x+1, thetaRes);
        for (int y = 0; y < phiRes; ++y) {
            const float t = (y + 0.5f) / phiRes,
                        phi = 2.f * M_PI * t;

            Vector v = SphericalDirection(sintheta, costheta, phi);
            gGrid[y*thetaRes + x] = Ashikhmin::computeGFactor(v, distribution);
            //gGrid[y*thetaRes + x] = std::max(1.f-2*s, 0.f); // TODO
        }
    }

    // TODO: test
    //mGFactorGrid = new MIPMap<float>(thetaRes, phiRes, gGrid.data(), false, 8.f, TEXTURE_CLAMP);

    mGFactorGrid.init(0, M_PI, thetaRes, 0, 2*M_PI, phiRes, gGrid);

    // TODO: implement in InterpolatedGrid?
    /*// Print g grid
    if (sPrintGGrid) {
        const int res = 50;
        fprintf(stderr, "g grid:\n");
        for (int x = 0; x < res; ++x) {
            const float s = (x + 0.5f) / res,
                        theta = M_PI * s,
                        costheta = cosf(theta),
                        sintheta = sinf(theta);
            for (int y = 0; y < res; ++y) {
                const float t = (y + 0.5f) / phiRes,
                            phi = 2.f * M_PI * t;

                Vector v = SphericalDirection(sintheta, costheta, phi);
                fprintf(stderr, " %.5ef", gFactor(v));
            }
            fprintf(stderr, "\n");
        }
    }*/
}

float
AshikhminCache::gFactor(const Vector &v) const
{
    //Assert(mGFactorGrid != NULL);
    const float theta = SphericalTheta(v),
                //s = theta * INV_PI,
                phi = SphericalPhi(v);
                //t = phi * INV_TWOPI;
    //return mGFactorGrid->Lookup(s, t);

    return mGFactorGrid.eval(theta, phi);
}

float
AshikhminCache::averageNH(void) const
{
    return mAvgNH;
}

const AshikhminCache&
AshikhminCache::get(const AshikhminDistribution &distribution)
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
                       const AshikhminDistribution &d)   // TODO: use reference
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
     R(reflectance), mDistribution(d), fresnel(f),
     mCache(AshikhminCache::get(d))
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
    float g_wi = gFactor(wi),
          g_wo = gFactor(wo);

    // TODO: show direct values:
    //       To show direct values without sampling, use a point light with white light,
    //       which the rendering equation becomes
    //         Ld += f * AbsDot(wi, n);  // See EstimateDirect() in integrator.cpp
    //       Note the wi and n is in world space; without bump mapping, in tangent space it should be wi.z
    //Spectrum ret;
    //float rgb[] = {SphericalTheta(wo)*INV_PI, SphericalPhi(wo)*INV_TWOPI, 0.f};
    //float rgb[] = {0.5f*(wo.x+1), 0.5f*(wo.y+1), 0.5f*(wo.z+1)};
    //ret = RGBSpectrum::FromRGB(rgb);
    //ret = g_wo;
    //ret = computeGFactor(wo, *mDistribution);
    /*Quaternion q(Vector(0, 0, 1), wo);
    Transform t = q.ToTransform();
    Vector v = t(Vector(0, 0, -1));
    float rgb[] = {0.5f*(v.x+1), 0.5f*(v.y+1), 0.5f*(v.z+1)};
    ret = RGBSpectrum::FromRGB(rgb);*/

    //ret /= wi.z;
    //return ret;

    // TODO: we need to make sure distribution->D(wh) actually returns a valid pdf; that is, it integrates to one over whole sphere
    return R * mDistribution.D(wh) * avgNH * F /
               (4.f * max(g_wi * g_wo, sSmallValue));
}

Spectrum
Ashikhmin::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const 
{
    if (mUseUniformSampling) {
        return BxDF::Sample_f(wo, wi, u1, u2, pdf);

    } else {
        if (wo.z < 0.f) {
            // Reverse side
            mDistribution.Sample_f(-wo, wi, u1, u2, pdf);
            *wi = -(*wi);
        } else {
            mDistribution.Sample_f(wo, wi, u1, u2, pdf);
        }

        if (!SameHemisphere(wo, *wi)) {
            return Spectrum(0.f);
        }
        return f(wo, *wi);
    }
}

float
Ashikhmin::Pdf(const Vector &wo, const Vector &wi) const 
{
    if (mUseUniformSampling) {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0.f;

    } else {
        if (!SameHemisphere(wo, wi)) {
            return 0.f;
        }

        if (wo.z < 0.f) {
            // Reverse side
            return mDistribution.Pdf(-wo, -wi);
        } else {
            return mDistribution.Pdf(wo, wi);
        }
    }
}

float
Ashikhmin::averageNH(void) const
{
    // TODO: test
    //return 1.f;
    //return computeAverageNH(mDistribution);
    return mCache.averageNH();
}

float
Ashikhmin::computeAverageNH(const AshikhminDistribution &distribution)
{
    // First dimension is phi, second theta
    const double xmin[2] = {0, 0},
                 xmax[2] = {2.f*M_PI, 0.5f*M_PI};
    double val = 0.f, error = 0.f;
#ifdef DEBUG
    int ret =
#endif
    hcubature(1, averageNHIntegrand, (void*)&distribution,
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
    const AshikhminDistribution *distribution =
                reinterpret_cast<const AshikhminDistribution*>(fdata);
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
    // TODO: test
    //return 1.f;
    //return computeGFactor(v, mDistribution);
    return mCache.gFactor(v);
}

float
Ashikhmin::computeGFactor(const Vector &v, const AshikhminDistribution &distribution)
{
    GFactorIntegrandData data;
    data.distribution = &distribution;
    data.v = v;

    // First dimension is phi, second theta
    // Note we now using method 2 in gFactorIntegrand(), so we integrate over the whole sphere
    const double xmin[2] = {0, 0},
                 xmax[2] = {2.f*M_PI, M_PI};
    double val = 0.f, error = 0.f;
#ifdef DEBUG
    int ret =
#endif
    hcubature(1, gFactorIntegrand, (void*)&data,
            2, xmin, xmax,
            0, sCubatureAbsError, sCubatureRelError,
            ERROR_INDIVIDUAL, &val, &error);

    Assert(ret == 0);

    //val = error;    // TODO: test

    // TODO: test
    //fprintf(stderr, "hcubature result (%d): %.8f +- %.8f\n",
    //        ret, val, error);

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

    // TODO: test
    // Method 1: treat the H as a vector in the reference frame determined by v
    //           so dot(v,H) = costheta, and we need to rotate v into the
    //           reference frame determined by N(0, 0, 1) before calling D()
    //           Note the integration domain is phi in [0,2*pi], theta in [0,pi/2]
    /*float VdotH = costheta;

    // Find H in the space defined by N
    Vector H_N = data->nToV(H);

    fval[0] = VdotH * data->distribution->D(H_N) * sintheta;*/

    // Method 2: alternatively, we can integrate over the whole sphere, and multiply
    //           an indicator function I(H) into the integrand:
    //             I(H) = 1 if dot(v,H) > 0
    //                    0 otherwise
    //           Note the H is now in the reference frame determined by N(0,0,1),
    //           so we need to really compute dot(v,H). The D() can now be evaluated
    //           directly against H. Notice the sintheta is unchanged, since it's only
    //           related to the parameterization, not the reference frame we used
    //           Note the integration domain is phi in [0,2*pi], theta in [0,pi]
    const float VdotH = Dot(data->v, H);
    if (VdotH > 0) {
        const float DofH = data->distribution->D(H);
        fval[0] = VdotH * DofH * sintheta;
    } else {
        fval[0] = 0.f;
    }

    // TODO: test
    //fprintf(stderr, "(%.2f, %.2f): %e\n", phi*INV_TWOPI, theta*INV_PI, fval[0]);

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
    // TODO: this is a wo direction where the first method in gFactorIntegrand() would fail
    //       (i.e., return 0)
    //const Vector VTest(-0.33629645195f, -0.62566874881f, 0.70387734241f);

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
