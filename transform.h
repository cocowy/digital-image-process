#ifndef TRANSFORM_H
#define TRANSFORM_H
#include <complex>
#include <QImage>
#include "fftw3.h"

enum ImageFilterType{
    IdealLowPass = 0,
    IdealHighPass,
    GaussainLowPass,
    GaussainHighPass,
    ButterworthLowPass,
    ButterworthHighPass,
};

QImage imageFFT2D(QImage src);
void imageFilterFFT2D(QImage src, int r, int option, QImage &originalSpectrumImage,
                      QImage &filteredSpectrumImage, QImage &dstImage);
void imageFilterFFT2D(fftwf_complex *y, int w, int h, float *filter,
                      QImage &filteredSpectrumImage, QImage &dstImage);
void generateFilter(int w, int h, int r, ImageFilterType type, float *filter);
void fftw2d(float *x, int w, int h, fftwf_complex *y);
void fftshift2D(fftwf_complex *src, int w, int h, fftwf_complex *dst);
void calcImageSpectrum(QImage src, QImage &dst);
void spectrum2QImage(fftwf_complex *s, int width, int height, QImage &dst);
void IdealLowPassFilter(QImage src, QImage &dst);
void IdealHighPassFilter(QImage src, QImage &dst);
void ButterworthLowPassFilter(QImage src, QImage &dst);
void ButterworthHighPassFilter(QImage src, QImage &dst);
void GaussainLowPassFilter(QImage src, QImage &dst);
void GaussainHighPassFilter(QImage src, QImage &dst);

#endif // TRANSFORM_H
