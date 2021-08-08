#include "imageprocess.h"

void paddingZeros(const float *r, const float *g, const float *b,
                   const int width, const int height,
                   const int half_pad_width, const int half_pad_height,
                   float *nr, float *ng, float *nb)
{
    int nw = width + 2*half_pad_width;
    int nh = height + 2*half_pad_height;

    memset(nr, 0, nh*nw*sizeof(float));
    memset(ng, 0, nh*nw*sizeof(float));
    memset(nb, 0, nh*nw*sizeof(float));

    for (int j=half_pad_height; j<nh-half_pad_height; j++) {
        memcpy(nr+j*nw+half_pad_width, r+(j-half_pad_height)*width, width*sizeof(float));
        memcpy(ng+j*nw+half_pad_width, g+(j-half_pad_height)*width, width*sizeof(float));
        memcpy(nb+j*nw+half_pad_width, b+(j-half_pad_height)*width, width*sizeof(float));
    }
}

void splitImageChannel(QImage &image, float *r, float *g, float *b)
{
    int width = image.width();
    int height = image.height();
    int count = 0;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = image.pixel(i, j);
            r[count] = qRed(pixel);
            g[count] = qGreen(pixel);
            b[count] = qBlue(pixel);
            count++;
        }
    }
}

void splitImageChannel(QImage &image, float *rgb)
{
    int width = image.width();
    int height = image.height();
    int count = 0;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = image.pixel(i, j);
            rgb[j*3*width+3*i] = (float)qRed(pixel);
            rgb[j*3*width+3*i+1] = (float)qGreen(pixel);
            rgb[j*3*width+3*i+2] = (float)qBlue(pixel);
            count++;
        }
    }
}

void splitImageChannel(QImage &image, uchar *r, uchar *g, uchar *b)
{
    int width = image.width();
    int height = image.height();
    int count = 0;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = image.pixel(i, j);
            r[count] = qRed(pixel);
            g[count] = qGreen(pixel);
            b[count] = qBlue(pixel);
            count++;
        }
    }
}

void splitImageChannel(QImage &image, uchar *rgb)
{
    int width = image.width();
    int height = image.height();
    int count = 0;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = image.pixel(i, j);
            rgb[j*3*width+3*i] = (uchar)qRed(pixel);
            rgb[j*3*width+3*i+1] = (uchar)qGreen(pixel);
            rgb[j*3*width+3*i+2] = (uchar)qBlue(pixel);
            count++;
        }
    }
}

void concatenateImageChannel(float *r, float *g, float *b, int w, int h, QImage &image)
{
    image = QImage(w, h, QImage::Format_RGB888);
    for (int j=0; j<h; j++) {
        for (int i=0; i<w; i++) {
            QRgb pixel = qRgb((int)r[j*w+i], (int)g[j*w+i], (int)b[j*w+i]);
            image.setPixel(i, j, pixel);
        }
    }
}

void concatenateImageChannel(float *rgb, int w, int h, QImage &image)
{
    image = QImage(w, h, QImage::Format_RGB888);
    for (int j=0; j<h; j++) {
        for (int i=0; i<w; i++) {
            int r = (int)rgb[j*3*w+3*i];
            int g = (int)rgb[j*3*w+3*i+1];
            int b = (int)rgb[j*3*w+3*i+2];
            QRgb pixel = qRgb(r, g, b);
            image.setPixel(i, j, pixel);
        }
    }
}

void concatenateImageChannel(uchar *r, uchar *g, uchar *b, int w, int h, QImage &image)
{
    image = QImage(w, h, QImage::Format_RGB888);
    for (int j=0; j<h; j++) {
        for (int i=0; i<w; i++) {
            QRgb pixel = qRgb((int)r, (int)g, (int)b);
            image.setPixel(i, j, pixel);
        }
    }
}

void concatenateImageChannel(uchar *rgb, int w, int h, QImage &image)
{
    image = QImage(w, h, QImage::Format_RGB888);
    for (int j=0; j<h; j++) {
        for (int i=0; i<w; i++) {
            int r = (int)rgb[j*3*w+3*i];
            int g = (int)rgb[j*3*w+3*i+1];
            int b = (int)rgb[j*3*w+3*i+2];
            QRgb pixel = qRgb(r, g, b);
            image.setPixel(i, j, pixel);
        }
    }
}

void rgb2ycrcb(uchar *r, uchar *g, uchar *b, int size, float *y, float *cr, float *cb)
{
    for (int i=0; i<size; i++)
    {
        y[i] = 0.256789 * r[i] + 0.504129 * g[i] + 0.097906 * b[i] + 16;
        cb[i] = -0.148223 * r[i] - 0.290992 * g[i] + 0.439215 * b[i] + 128;
        cr[i] = 0.439215 * r[i] - 0.367789 * g[i] - 0.071426 * b[i] + 128;
    }
}

void ycrcb2rgb(float *y, float *cr, float *cb, int size, uchar *r, uchar *g, uchar *b)
{
    for (int i=0; i<size; i++)
    {
        r[i] = 1.164383 * (y[i]-16) + 1.596027 * (cr[i]-128);
        g[i] = 1.164383 * (y[i]-16) - 0.391762 * (cb[i]-128)- 0.812969 * (cr[i]-128);
        b[i] = 1.164383 * (y[i]-16) + 2.017230 * (cb[i]-128);
    }
}

void qimage2ycrcb(QImage image, float *y, float *cr, float *cb)
{
    int pixel_num = image.width()*image.height();
    // obtain image channels
    uchar *channels = new uchar[pixel_num*3];
    uchar *r = channels;
    uchar *g = channels+pixel_num;
    uchar *b = channels+2*pixel_num;
    splitImageChannel(image, r, g, b);

    // rgb to ycrcb
    rgb2ycrcb(r, g, b, pixel_num, y, cr, cb);

    delete [] channels;
}

void ycrcb2qimage(float *y, float *cr, float *cb, int width, int height, QImage &image)
{
    int pixel_num = width*height;

    uchar *channels = new uchar[pixel_num*3];
    uchar *r = channels;
    uchar *g = channels+pixel_num;
    uchar *b = channels+2*pixel_num;

    // ycrcb to rgb
    ycrcb2rgb(y, cr, cb, pixel_num, r, g, b);

    // update image
    int count = 0;
    image = QImage(width, height, QImage::Format_RGB888);
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            int nr = r[count];
            int ng = g[count];
            int nb = b[count];
            count++;
            image.setPixel(i, j, qRgb(nr,ng,nb));
        }
    }
    delete [] channels;
}


QImage calculateHistogram(QImage &image, ImageChannel channel)
{
    // obtain gray image
    QImage grayImage = image.convertToFormat(QImage::Format_Grayscale8);
    int width = grayImage.width();
    int height = grayImage.height();
    uchar *r = new uchar[width*height];
    uchar *g = new uchar[width*height];
    uchar *b = new uchar[width*height];

    // seperate RGB channel
    splitImageChannel(image, r, g, b);

    QRgb hist_ior;

    uchar *bits = nullptr;
    switch (channel) {
        case ImageChannel::Y:
            hist_ior = qRgba(128, 128, 128, 255);
            bits = grayImage.bits();
            break;
        case ImageChannel::R:
            hist_ior = qRgba(255, 0, 0, 255);
            bits = r;
            break;
        case ImageChannel::G:
            hist_ior = qRgba(0, 255, 0, 255);
            bits = g;
            break;
        case ImageChannel::B:
            hist_ior = qRgba(0, 0, 255, 255);
            bits = b;
            break;
    }

    const int gray_level = 256;
    int hist[gray_level] = {0};

    // calculate histogram
    for (int i=0; i<width*height; i++)
    {
       int val = bits[i];
       hist[val]++;
    }

    // compress histogram into hist_image height
    int max_hist_val = hist[0];
    for (int i=1; i<gray_level; i++)
    {
        max_hist_val = hist[i] > max_hist_val ? hist[i] : max_hist_val;
    }

    int s_w = 2;
    int w = s_w*gray_level;
    int h = w;//width*height;
    float s_h = 0.8;

    for (int i=0; i<gray_level; i++)
    {
        int v = hist[i];
        hist[i] = int(h*1.0/max_hist_val * s_h * v);
    }

    QImage hist_image(w, h, QImage::Format_RGBA8888);
    QRgb value;
    for (int j=0; j<h; j++)
    {
        for (int i=0; i<w; i++)
        {
            if ( (hist[i/s_w] > 0) && (j >= h-hist[i/s_w]) )
                value = hist_ior;
            else
                value = qRgba(255, 255, 255, 255);
            hist_image.setPixel(i, j, value);
        }
    }
    delete [] r;
    delete [] g;
    delete [] b;

    return hist_image;
}

QImage calculateNegative(QImage &image, ImageChannel channel)
{
    // obtain gray image
    int width = image.width();
    int height = image.height();

    QImage newImage = image;
    int r, g, b;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = image.pixel(i, j);
            r = qRed(pixel);
            g = qGreen(pixel);
            b = qBlue(pixel);
            QRgb newPixel;
            switch (channel) {
                case ImageChannel::Y:
                    newPixel = qRgb(255-r, 255-g, 255-b);
                    break;
                case ImageChannel::R:
                    newPixel = qRgb(255-r, g, b);
                    break;
                case ImageChannel::G:
                    newPixel = qRgb(r, 255-g, b);
                    break;
                case ImageChannel::B:
                    newPixel = qRgb(r, g, 255-b);
                    break;
            }
            newImage.setPixel(i, j, newPixel);
          }
    }

    return newImage;
}

QImage convertToPseudoColor(QImage &image, ColorMap map)
{
    int width = image.width();
    int height = image.height();
    // obtain gray image
    QImage grayImage = image.convertToFormat(QImage::Format_Grayscale8);

    // seperate RGB channel
    QImage newImage(width, height, QImage::Format_RGB888);
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            QRgb pixel = grayImage.pixel(i, j);
            QRgb newPixel;
            int index, r, g, b;
            index = qRed(pixel); // gray scale
            r = qRed(pixel);
            g = qGreen(pixel);
            b = qBlue(pixel);
            switch (map) {
                case ColorMap::Jet:
                    r = jet_table[index*3];
                    g = jet_table[index*3+1];
                    b = jet_table[index*3+2];
                    break;
                case ColorMap::Parula:
                    r = parula_table[index*3];
                    g = parula_table[index*3+1];
                    b = parula_table[index*3+2];
                    break;
                case ColorMap::Hot:
                    r = hot_table[index*3];
                    g = hot_table[index*3+1];
                    b = hot_table[index*3+2];
                    break;
            }
            newPixel = qRgb(r, g, b);
            newImage.setPixel(i, j, newPixel);
          }
    }
    return newImage;
}

QImage equalizeHistogramProc1(QImage &image)
{
    int width = image.width();
    int height = image.height();
    int pixel_num = width*height;

    // obtain image channels
    uchar *channels = new uchar[pixel_num*3];
    uchar *r = channels;
    uchar *g = channels+pixel_num;
    uchar *b = channels+2*pixel_num;
    splitImageChannel(image, r, g, b);

    // rgb to ycrcb
    float *ycrcb = new float[pixel_num*3];
    float *y = ycrcb;
    float *cr = ycrcb+pixel_num;
    float *cb = ycrcb+2*pixel_num;
    rgb2ycrcb(r, g, b, pixel_num, y, cr, cb);

    // calculate hist/pdf
    int *hist = new int[pixel_num]; // hist/pdf
    const int gray_level = 256;
    float *gray_distribution = new float[gray_level];// cdf

    uchar *gray_equal = new uchar[gray_level]; // equalized gray

        // calculate pdf
        memset(hist, 0, pixel_num*sizeof(int));
        for (int i=0; i<pixel_num; i++)
        {
            int index = (int)y[i]; // gray scale
            hist[index]++;
        }

        // calculate cdf
        memset(gray_distribution, 0, gray_level*sizeof(float));
        gray_distribution[0] = hist[0]*1.0f/pixel_num;
        for (int i = 1; i < gray_level; i++)
        {
            gray_distribution[i] = gray_distribution[i-1] + hist[i]*1.0f/pixel_num;
        }

        // recalculate equalized gray
        memset(gray_equal, 0, gray_level*sizeof(uchar));
        for (int i = 0; i < gray_level; i++)
        {
            gray_equal[i] = (uchar)(255 * gray_distribution[i] + 0.5);
        }

        // new gray channel
        for (int i=0; i<pixel_num; i++)
        {
            int index = (int)y[i]; // gray scale
            y[i] = gray_equal[index];
        }

        // ycrcb to rgb
    ycrcb2rgb(y, cr, cb, pixel_num, r, g, b);

    // update image
    int count = 0;
    QImage newImage = image;//grayImage;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            int nr = r[count];
            int ng = g[count];
            int nb = b[count];
            count++;
            newImage.setPixel(i, j, qRgb(nr,ng,nb));
        }
    }
    delete [] gray_equal;
    delete [] gray_distribution;
    delete [] channels;
    delete [] ycrcb;

    return newImage;
}

QImage equalizeHistogramProc(QImage &image)
{
    int width = image.width();
    int height = image.height();
    int pixel_num = width*height;

    // obtain gray image
    uchar *channels = new uchar[width*height*3];
    uchar *r = channels;
    uchar *g = channels+width*height;
    uchar *b = channels+2*width*height;
    splitImageChannel(image, r, g, b);

    uchar *c[4] = {r, g, b, 0};

    // calculate hist/pdf
    int *hist = new int[pixel_num]; // hist/pdf
    const int gray_level = 256;
    float *gray_distribution = new float[gray_level];// cdf

    uchar *gray_equal = new uchar[gray_level]; // equalized gray
    for (uchar **p=c; (*p) != 0; p++)
    {
        // calculate pdf
        memset(hist, 0, pixel_num*sizeof(int));
        for (int i=0; i<pixel_num; i++)
        {
            int index = (*p)[i]; // gray scale
            hist[index]++;
        }

        // calculate cdf
        memset(gray_distribution, 0, gray_level*sizeof(float));
        gray_distribution[0] = hist[0]*1.0f/pixel_num;
        for (int i = 1; i < gray_level; i++)
        {
            gray_distribution[i] = gray_distribution[i-1] + hist[i]*1.0f/pixel_num;
        }

        // recalculate equalized gray
        memset(gray_equal, 0, gray_level*sizeof(uchar));
        for (int i = 0; i < gray_level; i++)
        {
            gray_equal[i] = (uchar)(255 * gray_distribution[i] + 0.5);
        }

        // new gray channel
        for (int i=0; i<pixel_num; i++)
        {
            int index = (*p)[i]; // gray scale
            (*p)[i] = gray_equal[index];
        }
    }

    // update image
    int count = 0;
    QImage newImage = image;//grayImage;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            int nr = r[count];
            int ng = g[count];
            int nb = b[count];
            count++;
            newImage.setPixel(i, j, qRgb(nr,ng,nb));
        }
    }
    delete [] gray_equal;
    delete [] gray_distribution;
    delete [] channels;

    return newImage;
}

void filter(float *src, int image_width, int image_height, int coefficient,
            float *kernel, int kernel_width, int kernel_height, float* dst)
{
    int i, j, m, n;
    float val;

    memset(dst, 0, image_width*image_height);

    for (j=kernel_height/2; j<image_height-kernel_height/2; j++)
    {
        for (i=kernel_width/2; i<image_width-kernel_width/2; i++)
        {
            val = 0;
            for (n=-kernel_height/2; n<=kernel_height/2; n++)
            {
                for (m=-kernel_width/2; m<=kernel_width/2; m++)
                {
                    val += src[(j-n)*image_width+(i-m)] *
                            kernel[(n+kernel_height/2)*kernel_width+m+kernel_width/2];
                }
            }
            val = val/coefficient;
            val = val>255 ? 255 : val;
            dst[j*image_width+i] = (float) (val<0? 0:val);
        }
    }
}

void edgefilter(float *src, int image_width, int image_height, float *kernel1, float *kernel2, int kernel_width,
                int kernel_height, float* dst)
{
    int i, j, m, n;
    float val;
    val=0;
    float* tmp= src;
    memset(dst, 0, image_width*image_height);

    for (j=kernel_height/2; j<image_height-kernel_height/2; j++)
    {
        for (i=kernel_width/2; i<image_width-kernel_width/2; i++)
        {
            val = 0;
            for (n=-kernel_height/2; n<kernel_height/2; n++)
            {
                for (m=-kernel_width/2; m<=kernel_width/2; m++)
                {
                    val += src[(j-n)*image_width+(i-m)] *
                            kernel1[(n+kernel_height/2)*kernel_width+m+kernel_width/2];
               }
            }
            tmp[j*image_width+i] = abs(val);
        }
    }

    for (j=kernel_height/2; j<image_height-kernel_height/2; j++)
    {
        for (i=kernel_width/2; i<image_width-kernel_width/2; i++)
        {
            val = 0;
            for (n=-kernel_height/2; n<kernel_height/2; n++)
            {
                for (m=-kernel_width/2; m<=kernel_width/2; m++)
                {
                    val += tmp[(j-n)*image_width+(i-m)] * kernel2[(n+kernel_height/2)*kernel_width+m+kernel_width/2];
                }
            }
            dst[j*image_width+i] = abs(val);
        }
    }
}
void logfilter(float *src, int image_width, int image_height, float *kernel, int kernel_width, int kernel_height, float* dst)
{
    int i, j, m, n;
    float val;
    val=0;
    memset(dst, 0, image_width*image_height);


    for (j=kernel_height/2; j<image_height-kernel_height/2; j++)
    {
        for (i=kernel_width/2; i<image_width-kernel_width/2; i++)
        {
            val = 0;
            for (n=-kernel_height/2; n<kernel_height/2; n++)
            {
                for (m=-kernel_width/2; m<=kernel_width/2; m++)
                {
                    val += src[(j-n)*image_width+(i-m)] * kernel[(n+kernel_height/2)*kernel_width+m+kernel_width/2];
                }
            }
            dst[j*image_width+i] = abs(val);
        }
    }
}

///edge filtering
void EdgeFilterRobert(QImage src, QImage &dst)
{
    int width = src.width();
    int height = src.height();
    float kernel1[] = {1,0,0,-1};
    float kernel2[] = {0,1,-1,0};

    int kernel_width = 2;
    int kernel_height = 2;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(src, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    filter(R, width, height, 1, kernel1, kernel_width, kernel_height, dst1);
    filter(G, width, height, 1, kernel1, kernel_width, kernel_height, dst2);
    filter(B, width, height, 1, kernel1, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, dst);
}
void EdgeFilterSobel(QImage src, QImage &dst)
{
    int width = src.width();
    int height = src.height();
    float kernel1[] = {-1,-2,-1,0,0,0,1,2,1};
    float kernel2[] = {-1,0,1,-2,0,2,-1,0,1};

    int kernel_width = 3;
    int kernel_height = 3;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(src, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    edgefilter(R, width, height, kernel1, kernel2, kernel_width, kernel_height, dst1);
    edgefilter(G, width, height, kernel1, kernel2, kernel_width, kernel_height, dst2);
    edgefilter(B, width, height, kernel1, kernel2, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, dst);
}
void EdgeFilterLOG(QImage src, QImage &dst)
{
    QImage tmp = src;
    int width = src.width();
    int height = src.height();
    float kernel1[] = {1,1,1,1,-8,1,1,1,1};

    int kernel_width = 3;
    int kernel_height = 3;

    GaussainLowPassFilter(src, tmp);
    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(tmp, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    logfilter(R, width, height, kernel1, kernel_width, kernel_height, dst1);
    logfilter(G, width, height, kernel1, kernel_width, kernel_height, dst2);
    logfilter(B, width, height, kernel1, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, dst);
}

void EdgeFilterPerwitt(QImage src, QImage &dst)
{
    int width = src.width();
    int height = src.height();
    float kernel1[] = {-1,-1,-1,0,0,0,1,1,1};
    float kernel2[] = {-1,0,1,-1,0,1,-1,0,1};

    int kernel_width = 3;
    int kernel_height = 3;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(src, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    edgefilter(R, width, height,  kernel1, kernel2, kernel_width, kernel_height, dst1);
    edgefilter(G, width, height,  kernel1, kernel2, kernel_width, kernel_height, dst2);
    edgefilter(B, width, height,  kernel1, kernel2, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, dst);
}


////spatial filtering
QImage AverageFilter(QImage &image) //3*3 average filter
{
    QImage newImage = image;//grayImage;
    int width = image.width();
    int height = image.height();
    int len = 9;
    float* kernel = (float*)malloc(len*sizeof(len));
    std::fill_n(kernel, len, 1);
    int kernel_width = 3;
    int kernel_height = 3;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(image, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    filter(R, width, height, 9, kernel, kernel_width, kernel_height, dst1);
    filter(G, width, height, 9, kernel, kernel_width, kernel_height, dst2);
    filter(B, width, height, 9, kernel, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, newImage);
    return newImage;
}

QImage GaussianFilter(QImage &image) //3*3 Gaussian filter
{
    QImage newImage = image;//grayImage;
    int width = image.width();
    int height = image.height();
    //int len = 9;
    //float* kernel = (float*)malloc(len*sizeof(len));
    float kernel[] = {1,2,1,2,4,2,1,1,1};
    int kernel_width = 3;
    int kernel_height = 3;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(image, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    filter(R, width, height, 16, kernel, kernel_width, kernel_height, dst1);
    filter(G, width, height, 16, kernel, kernel_width, kernel_height, dst2);
    filter(B, width, height, 16, kernel, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, newImage);
    return newImage;
}


float Sort(float a1[],int n)
{
    int k;
    int j;
    int c;
    for (int i = 0;i<n;i++)
    {
        k = i;
        for(j = k +1;j<n;j++)
        {
        if (a1[i]>a1[j])
        {
            c = a1[j];
            a1[j] = a1[i];
            a1[i] = c;
        }
        }
    }
    return a1[n/2];
}


void median_filter(float *src, int image_width, int image_height, int coefficient,
            float *kernel, int kernel_width, int kernel_height, float* dst)
{
    int i, j, m, n;
    float* val=new float [image_height*image_width];

    memset(dst, 0, image_width*image_height);

    for (j=kernel_height/2; j<image_height-kernel_height/2; j++)
    {
        for (i=kernel_width/2; i<image_width-kernel_width/2; i++)
        {
            int cnt = 0;
            for (n=-kernel_height/2; n<=kernel_height/2; n++)
            {
                for (m=-kernel_width/2; m<=kernel_width/2; m++)
                {
                    val[cnt] = src[(j-n)*image_width+(i-m)] *
                            kernel[(n+kernel_height/2)*kernel_width+m+kernel_width/2];
                    cnt++;

                }
            }
            int res = Sort(val, 9);
            res = res/coefficient;
            res = res>255 ? 255 : res;
            dst[j*image_width+i] = (float) (res<0? 0:res);
        }
    }
}

QImage MedianFilter(QImage &image)//3*3 median filter
{
    QImage newImage = image;//grayImage;
    int width = image.width();
    int height = image.height();
    //int len = 9;
    //float* kernel = (float*)malloc(len*sizeof(len));
    float kernel[] = {1,1,1,1,1,1,1,1,1};
    int kernel_width = 3;
    int kernel_height = 3;

    float* R = new float[width*height];
    float* G = new float[width*height];
    float* B = new float[width*height];
    splitImageChannel(image, R, G, B);

    float* dst1 = new float[width*height];
    float* dst2 = new float[width*height];
    float* dst3 = new float[width*height];

    median_filter(R, width, height, 1, kernel, kernel_width, kernel_height, dst1);
    median_filter(G, width, height, 1, kernel, kernel_width, kernel_height, dst2);
    median_filter(B, width, height, 1, kernel, kernel_width, kernel_height, dst3);
    concatenateImageChannel(dst1, dst2, dst3, width, height, newImage);
    return newImage;
}

QImage mat2qim(const cv::Mat& mat)
{
    // 8-bits unsigned, NO. OF CHANNELS = 1
    if(mat.type() == CV_8UC1)
    {
        QImage image(mat.cols, mat.rows, QImage::Format_Indexed8);
        // Set the color table (used to translate colour indexes to qRgb values)
        image.setColorCount(256);
        for(int i = 0; i < 256; i++)
        {
            image.setColor(i, qRgb(i, i, i));
        }
        // Copy input Mat
        uchar *pSrc = mat.data;
        for(int row = 0; row < mat.rows; row ++)
        {
            uchar *pDest = image.scanLine(row);
            memcpy(pDest, pSrc, mat.cols);
            pSrc += mat.step;
        }
        return image;
    }
    // 8-bits unsigned, NO. OF CHANNELS = 3
    else if(mat.type() == CV_8UC3)
    {
        // Copy input Mat
        const uchar *pSrc = (const uchar*)mat.data;
        // Create QImage with same dimensions as input Mat
        QImage image(pSrc, mat.cols, mat.rows, mat.step, QImage::Format_RGB888);
        return image.rgbSwapped();
    }
    else if(mat.type() == CV_8UC4)
    {
        qDebug() << "CV_8UC4";
        // Copy input Mat
        const uchar *pSrc = (const uchar*)mat.data;
        // Create QImage with same dimensions as input Mat
        QImage image(pSrc, mat.cols, mat.rows, mat.step, QImage::Format_ARGB32);
        return image.copy();
    }
    else
    {
        qDebug() << "ERROR: Mat could not be converted to QImage.";
        return QImage();
    }
}

void qimage_to_mat(const QImage& image, cv::OutputArray out) {

    switch(image.format()) {
        case QImage::Format_Invalid:
        {
            cv::Mat empty;
            empty.copyTo(out);
            break;
        }
        case QImage::Format_RGB32:
        {
            cv::Mat view(image.height(),image.width(),CV_8UC4,(void *)image.constBits(),image.bytesPerLine());
            view.copyTo(out);
            break;
        }
        case QImage::Format_RGB888:
        {
            cv::Mat view(image.height(),image.width(),CV_8UC3,(void *)image.constBits(),image.bytesPerLine());
            cvtColor(view, out, cv::COLOR_RGB2BGR);
            break;
        }
        default:
        {
            QImage conv = image.convertToFormat(QImage::Format_ARGB32);
            cv::Mat view(conv.height(),conv.width(),CV_8UC4,(void *)conv.constBits(),conv.bytesPerLine());
            view.copyTo(out);
            break;
        }
    }
}


void SIFTMatching(QImage src, QImage &dst)
{
    int numFeatures = 100;
    cv::Ptr<cv::Feature2D> sift = cv::SiftFeatureDetector::create(numFeatures);
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat srcImage;
    qimage_to_mat(src, srcImage);
    sift->detect(srcImage,keyPoints);

    cv::Mat keypoint_img;
    cv::drawKeypoints(srcImage, keyPoints, keypoint_img, cv::Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    dst=mat2qim(keypoint_img);
}

void FastMatching(QImage src, QImage &dst)
{

    int numFeatures = 40;
    cv::Ptr<cv::Feature2D> fast = cv::FastFeatureDetector::create(numFeatures);
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat srcImage;
    qimage_to_mat(src, srcImage);
    fast->detect(srcImage,keyPoints, Mat());

    Mat keypoint_img;
    cv::drawKeypoints(srcImage, keyPoints, keypoint_img, Scalar::all(-1), DrawMatchesFlags::DEFAULT);
    dst=mat2qim(keypoint_img);
}

void ORBMatching(QImage src, QImage &dst)
{
    int numFeatures = 100;
    cv::Ptr<cv::Feature2D> orb = cv::ORB::create(numFeatures);
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat srcImage;
    qimage_to_mat(src, srcImage);
    orb->detect(srcImage,keyPoints, Mat());

    Mat keypoint_img;
    drawKeypoints(srcImage, keyPoints, keypoint_img, Scalar::all(-1), DrawMatchesFlags::DEFAULT);
    dst=mat2qim(keypoint_img);
}


void AddNoiseSaltPepper(QImage src, QImage &dst)
{
    cv::Mat srcImage;
    qimage_to_mat(src, srcImage);

    int n=5000;
    int c =srcImage.cols;
    int t =srcImage.rows;

    for(int k=0;k<n;k++)
    {
        int j = std::rand()%t;
        int i = std::rand()%c;
        if(srcImage.channels() == 1){
            srcImage.at<uchar>(j,i) = 255;
        }
        else{
            srcImage.at<cv::Vec3b>(j,i)[0] = 255;
            srcImage.at<cv::Vec3b>(j,i)[1] = 255;
            srcImage.at<cv::Vec3b>(j,i)[2] = 255;
        }
    }
    for(int k=0;k<n;k++)
    {
        int i = std::rand()%c;
        int j = std::rand()%t;
        if(srcImage.channels() == 1){
            srcImage.at<uchar>(j,i) = 0;
        }
        else{
            srcImage.at<cv::Vec3b>(j,i)[0] = 0;
            srcImage.at<cv::Vec3b>(j,i)[1] = 0;
            srcImage.at<cv::Vec3b>(j,i)[2] = 0;
        }
    }
    dst=mat2qim(srcImage);
}
void AddNoiseGuassian(QImage src, QImage &dst)
{
    cv::Mat srcImage;
    qimage_to_mat(src, srcImage);

    cv::Mat noise(srcImage.size(),srcImage.type());
    float m = (10,10,10);
    float sigma = (10,50,50);
    cv::randn(noise, m, sigma); //mean and variance
    srcImage += noise;
    dst=mat2qim(srcImage);
}

void adaptiveContrastEnhancement(QImage &src_image, float *rgb, float *rgb_ii, float *rgb_ii_power, int max_window_size,
                                 int half_window_size, float alpha, float max_cg, QImage &dst_image)
{
    dst_image = src_image;
    int image_width = src_image.width();
    int image_height = src_image.height();
    int pixel_num = image_width*image_height;

    int max_image_width = src_image.width() + 2*max_window_size;
    int max_image_height = src_image.height() + 2*max_window_size;
    int max_pixel_num = max_image_width*max_image_height;
    int max_kernel_height = 2*max_window_size+1;
    int max_kernel_width = 2*max_window_size+1;

    int i=0, j=0;
    int kernel_height = 2*half_window_size+1;
    int kernel_width = 2*half_window_size+1;
    int kernel_size = kernel_height*kernel_width;
    float image_mean=0, image_std=0;
    for (int c=0; c<3; c++)
    {
        // image mean
        image_mean = box_integral(rgb_ii+c*max_pixel_num, max_image_width, max_image_height,
                               max_window_size, max_window_size + image_width-1,
                               max_window_size, max_window_size + image_height-1);
        image_mean /= pixel_num;

        // image std
        image_std = box_integral(rgb_ii_power+c*max_pixel_num, max_image_width, max_image_height,
                               max_window_size, max_window_size + image_width-1,
                               max_window_size, max_window_size + image_height-1);
        image_std /= pixel_num;
        image_std -= image_mean*image_mean;

        image_std = sqrtf(image_std);

        // local area mean and std
        for (j=max_kernel_height/2; j<max_image_height-max_kernel_height/2; j++)
        {
            for (i=max_kernel_width/2; i<max_image_width-max_kernel_width/2; i++)
            {
                // mean
                float mean = box_integral(rgb_ii+c*max_pixel_num, max_image_width, max_image_height,
                                       i-kernel_width/2, i+kernel_width/2,
                                       j-kernel_height/2, j+kernel_height/2);
                mean /= kernel_size;

                // std
                float std= box_integral(rgb_ii_power+c*max_pixel_num, max_image_width, max_image_height,
                                       i-kernel_width/2, i+kernel_width/2,
                                       j-kernel_height/2, j+kernel_height/2);
                std = std/kernel_size - mean*mean;
                std = sqrtf(std);

                // constrast gain
                float cg = alpha*image_std/std;
                if (cg>max_cg) cg = max_cg;

                float dst_val = mean + cg * (rgb[c*max_pixel_num + j*max_image_width+i] - mean);
                if (dst_val > 255) dst_val = 255;
                if (dst_val < 0) dst_val = 0;
                QRgb temp = dst_image.pixel(i-max_kernel_width/2, j-max_kernel_height/2);
                int temp_r = qRed(temp);
                int temp_g = qGreen(temp);
                int temp_b = qBlue(temp);
                if (c==0) temp_r = dst_val;
                if (c==1) temp_g = dst_val;
                if (c==2) temp_b = dst_val;
                dst_image.setPixel(i-max_kernel_width/2, j-max_kernel_height/2, qRgb(temp_r, temp_g, temp_b));
            }
        }
    }
}

void calculate_integral_image(float *image, int width, int height, float *integral_image)
{
    int i, j;

    // first row only
    float rs = 0;
    for(j=0; j<width; j++)
    {
        rs += image[j];
        integral_image[j] = rs;
    }
    for(i=1; i<height; ++i)
    {
        rs = 0;
        for(j=0; j<width; ++j)
        {
            rs += image[i*width+j];
            integral_image[i*width+j] = rs + integral_image[(i-1)*width+j];
        }
    }
}

void calculate_integral_image_power(float *image, int width, int height, float *integral_image)
{
    int i, j;

    // first row only
    float rs = 0;
    for(j=0; j<width; j++)
    {
        rs += image[j]*image[j];
        integral_image[j] = rs;
    }
    for(i=1; i<height; ++i)
    {
        rs = 0;
        for(j=0; j<width; ++j)
        {
            rs += image[i*width+j]*image[i*width+j];
            integral_image[i*width+j] = rs + integral_image[(i-1)*width+j];
        }
    }
}
__inline float box_integral(float *integral_image, int width, int height, int c1, int c2, int r1, int r2)
{
    float a, b, c, d;

    a = (c1-1<0 || r1-1<0)	? 0 : integral_image[(r1-1) * width + (c1-1)];
    b = r1-1<0				? 0 : integral_image[(r1-1) * width + c2];
    c = c1-1<0				? 0 : integral_image[r2 * width + (c1-1)];
    d = integral_image[r2 * width + c2];

    return a - b - c + d;
}
