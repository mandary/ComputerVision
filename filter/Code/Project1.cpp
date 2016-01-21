#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>

/***********************************************************************
  This is the only file you need to change for your assignment.  The
  other files control the UI (in case you want to make changes.)
************************************************************************/


// The first four functions provide example code to help get you started

// Convert an image to grey-scale
void MainWindow::BlackWhiteImage(QImage *image)
{
    int r, c;
    QRgb pixel;

    for(r=0;r<image->height();r++)
    {
        for(c=0;c<image->width();c++)
        {
            pixel = image->pixel(c, r);
            double red = (double) qRed(pixel);
            double green = (double) qGreen(pixel);
            double blue = (double) qBlue(pixel);

            // Compute intensity from colors - these are common weights
            double intensity = 0.3*red + 0.6*green + 0.1*blue;

            image->setPixel(c, r, qRgb( (int) intensity, (int) intensity, (int) intensity));
        }
    }
}

// Add random noise to the image
void MainWindow::AddNoise(QImage *image, double mag, bool colorNoise)
{
    int r, c;
    QRgb pixel;
    int noiseMag = mag;
    noiseMag *= 2;

    for(r=0;r<image->height();r++)
    {
        for(c=0;c<image->width();c++)
        {
            pixel = image->pixel(c, r);
            int red = qRed(pixel);
            int green = qGreen(pixel);
            int blue = qBlue(pixel);

            // If colorNoise, add color independently to each channel
            if(colorNoise)
            {
                red += rand()%noiseMag - noiseMag/2;
                green += rand()%noiseMag - noiseMag/2;
                blue += rand()%noiseMag - noiseMag/2;
            }
            // otherwise add the same amount of noise to each channel
            else
            {
                int noise = rand()%noiseMag - noiseMag/2;

                red += noise;
                green += noise;
                blue += noise;
            }

            // Make sure we don't over or under saturate
            red = min(255, max(0, red));
            green = min(255, max(0, green));
            blue = min(255, max(0, blue));

            image->setPixel(c, r, qRgb( red, green, blue));
        }
    }
}

// Here is an example of blurring an image using a mean or box filter with the specified radius.
// This could be implemented using separable filters to make it much more efficient, but it is not.
void MainWindow::MeanBlurImage(QImage *image, int radius)
{
    if(radius == 0)
        return;

    int r, c, rd, cd, i;
    QRgb pixel;

    // This is the size of the kernel
    int size = 2*radius + 1;

    // Create a buffer image so we're not reading and writing to the same image during filtering.
    QImage buffer;
    int w = image->width();
    int h = image->height();

    // This creates an image of size (w + 2*radius, h + 2*radius) with black borders.
    // This could be improved by filling the pixels using a different padding technique (reflected, fixed, etc.)
    buffer = image->copy(-radius, -radius, w + 2*radius, h + 2*radius);

    // Compute kernel to convolve with the image.
    double *kernel = new double [size*size];

    for(i=0;i<size*size;i++)
    {
        kernel[i] = 1.0;
    }

    // Make sure kernel sums to 1
    double denom = 0.000001;
    for(i=0;i<size*size;i++)
        denom += kernel[i];
    for(i=0;i<size*size;i++)
        kernel[i] /= denom;

    // For each pixel in the image...
    for(r=0;r<h;r++)
    {
        for(c=0;c<w;c++)
        {
            double rgb[3];

            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;

            // Convolve the kernel at each pixel
            for(rd=-radius;rd<=radius;rd++)
                for(cd=-radius;cd<=radius;cd++)
                {
                     // Get the pixel value
                     pixel = buffer.pixel(c + cd + radius, r + rd + radius);

                     // Get the value of the kernel
                     double weight = kernel[(rd + radius)*size + cd + radius];

                     rgb[0] += weight*(double) qRed(pixel);
                     rgb[1] += weight*(double) qGreen(pixel);
                     rgb[2] += weight*(double) qBlue(pixel);
                }

            // Store mean pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

    // Clean up.
    delete [] kernel;
}

// Downsample the image by 1/2
void MainWindow::HalfImage(QImage &image)
{
    QImage buffer;
    int w = image.width();
    int h = image.height();
    int r, c;

    buffer = image.copy();

    // Reduce the image size.
    image = QImage(w/2, h/2, QImage::Format_RGB32);

    // Copy every other pixel
    for(r=0;r<h/2;r++)
        for(c=0;c<w/2;c++)
        {
             image.setPixel(c, r, buffer.pixel(c*2, r*2));
        }
}



void MainWindow::GaussianBlurImage(QImage *image, double sigma)
{
    if(sigma == 0) // nothing if sigma is 0
        return;

    int r, c, rd, cd, i;
    QRgb pixel;

    // This is the size of the kernel
    int isig = (int)sigma;
    int size = 2 * isig + 1;

    // Create a buffer image so we're not reading and writing to the same image during filtering.
    QImage buffer;
    int w = image->width();
    int h = image->height();

    // This creates an image of size (w + 2*sigma, h + 2*sigma) with black borders.
    // This could be improved by filling the pixels using a different padding technique (reflected, fixed, etc.)
    buffer = image->copy(-isig, -isig, w + 2 * isig, h + 2 * isig);

    // Compute kernel to convolve with the image.
    double *kernel = new double [size * size];
    double sum = 0.000001;
    double mean = size / 2.0;
    for(i = 0; i < size * size; i++) {
        double twoSigSq = 2 * sigma * sigma;
        kernel[i] = exp(-(pow((i / size - mean), 2.0) + pow((i % size - mean), 2.0)) / twoSigSq)
                         / (M_PI * twoSigSq);
        sum += kernel[i];
    }

    // Normalize kernel weights
    for(i = 0; i < size * size; i++)
        kernel[i] /= sum;

    // Convolve kernel around image
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double rgb[3];

            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;

            // Convolve the kernel at each pixel
            for(rd = -isig; rd <= isig; rd++) {
                for(cd = -isig; cd <= isig; cd++) {
                     // Get the pixel value
                     pixel = buffer.pixel(c + cd + isig, r + rd + isig);

                     // Get the value of the kernel
                     double weight = kernel[(rd + isig) * size + cd + isig];

                     rgb[0] += weight*(double) qRed(pixel);
                     rgb[1] += weight*(double) qGreen(pixel);
                     rgb[2] += weight*(double) qBlue(pixel);
                }
            }

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

    // Clean up.
    delete [] kernel;
}

void MainWindow::SeparableGaussianBlurImage(QImage *image, double sigma)
{
    if(sigma == 0)
        return;

    int r, c, rd, cd, i;
    QRgb pixel;

    // This is the size of the kernel
    int isig = (int)sigma;
    int size = 2 * isig + 1;


    // Create a buffer image so we're not reading and writing to the same image during filtering.
    QImage buffer;
    int w = image->width();
    int h = image->height();

    // This creates an image of size (w + 2*sigma, h + 2*sigma) with black borders.
    // This could be improved by filling the pixels using a different padding technique (reflected, fixed, etc.)
    buffer = image->copy(-isig, -isig, w + 2 * isig, h + 2 * isig);

    // Compute kernel to convolve with the image.
    double* kernel = new double [size];
    double mean = size / 2.0;
    double sum = 0.000001;
    for(i = 0; i < size; i++) {
        kernel[i] = exp(-0.5 * pow((i - mean) / sigma, 2.0)) / (sqrt(2 * M_PI) * sigma);
        sum += kernel[i];
    }

    // Normalize kernel weights
    for(i = 0; i < size; i++)
        kernel[i] /= sum;

    // Convolve kernel around image
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double rgb[3];

            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;

            // Convolve the kernel at each pixel
            for(rd = -isig; rd <= isig; rd++) {
            
                pixel = buffer.pixel(c + isig, r + rd + isig);

                // Get the value of the kernel
                double weight = kernel[rd + isig];

                rgb[0] += weight*(double) qRed(pixel);
                rgb[1] += weight*(double) qGreen(pixel);
                rgb[2] += weight*(double) qBlue(pixel);
            }

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

    // make a new copy from updated image
    buffer = image->copy(-isig, -isig, w + 2 * isig, h + 2 * isig);

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double rgb[3];

            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;

            // Convolve the kernel at each pixel
            for(cd = -isig; cd <= isig; cd++) {
            
                pixel = buffer.pixel(c + cd + isig, r + isig);

                // Get the value of the kernel
                double weight = kernel[cd + isig];

                rgb[0] += weight*(double) qRed(pixel);
                rgb[1] += weight*(double) qGreen(pixel);
                rgb[2] += weight*(double) qBlue(pixel);
            }

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

    // Clean up.
    delete [] kernel;

}

void MainWindow::FirstDerivImage(QImage *image, double sigma)
{
    // For image derivatives really is (pixel - next pixel) / 1
    int r, c;

    QRgb pixel;

    QImage buffer;
    int w = image->width();
    int h = image->height();

    // Copy image
    buffer = image->copy(-1, -1, w + 2, h + 2);


    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            

            // Calibrate for negative value
            double rgb[3] = {128, 128, 128};
            

            pixel = buffer.pixel(c, r);

            rgb[0] -= (double) qRed(pixel);
            rgb[1] -= (double) qGreen(pixel);
            rgb[2] -= (double) qBlue(pixel);

            pixel = buffer.pixel(c + 2, r);

            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            // Adjust the rgb value
            rgb[0] = min(255.0, max(0.0, rgb[0]));
            rgb[1] = min(255.0, max(0.0, rgb[1]));
            rgb[2] = min(255.0, max(0.0, rgb[2]));

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));

        }
    }

    // Gaussian Blur
    GaussianBlurImage(image, sigma);

}

void MainWindow::SecondDerivImage(QImage *image, double sigma)
{
    // Again, derivative is really the difference
    // This time is all neighbors - center
    int r, c;

    QRgb pixel;

    QImage buffer;
    int w = image->width();
    int h = image->height();

    // Copy image
    buffer = image->copy(-1, -1, w + 2, h + 2);

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {

            // Calibrate for negative value
            double rgb[3] = {128, 128, 128};

            // Subtract itself
            pixel = buffer.pixel(c + 1, r + 1);
            rgb[0] -= 4 * (double) qRed(pixel);
            rgb[1] -= 4 * (double) qGreen(pixel);
            rgb[2] -= 4 * (double) qBlue(pixel);


            // Add neighbors
            pixel = buffer.pixel(c, r + 1);
            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            pixel = buffer.pixel(c + 2, r + 1);
            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            pixel = buffer.pixel(c + 1, r);
            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            pixel = buffer.pixel(c + 1, r + 2);
            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            // Adjust the rgb value
            rgb[0] = min(255.0, max(0.0, rgb[0]));
            rgb[1] = min(255.0, max(0.0, rgb[1]));
            rgb[2] = min(255.0, max(0.0, rgb[2]));

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));

        }
    }

    // Gaussian Blur
    GaussianBlurImage(image, sigma);
}

void MainWindow::SharpenImage(QImage *image, double sigma, double alpha)
{
    // sharpen = original - alpha * second derivative

    int r, c;

    QRgb pixel;

    QImage buffer;
    buffer = image->copy();

    int w = image->width();
    int h = image->height();

    SecondDerivImage(&buffer, sigma);

    for(r = 1; r < h - 1; r++) {
        for(c = 1; c < w - 1; c++) {

            double rgb[3] = {0, 0, 0};

            pixel = image->pixel(c, r);
            rgb[0] += (double) qRed(pixel);
            rgb[1] += (double) qGreen(pixel);
            rgb[2] += (double) qBlue(pixel);

            pixel = buffer.pixel(c, r);
            rgb[0] -= alpha * ((double) qRed(pixel) - 128.0);
            rgb[1] -= alpha * ((double) qGreen(pixel) - 128.0);
            rgb[2] -= alpha * ((double) qBlue(pixel) - 128.0);

            // Adjust the rgb value
            rgb[0] = min(255.0, max(0.0, rgb[0]));
            rgb[1] = min(255.0, max(0.0, rgb[1]));
            rgb[2] = min(255.0, max(0.0, rgb[2]));

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

}

void MainWindow::BilateralImage(QImage *image, double sigmaS, double sigmaI)
{
    // Add your code here.  Should be similar to GaussianBlurImage.
}

void MainWindow::SobelImage(QImage *image)
{
    // Add your code here.

    /***********************************************************************
      When displaying the orientation image I
      recommend the following:

    double mag; // magnitude of the gradient
    double orien; // orientation of the gradient

    double red = (sin(orien) + 1.0)/2.0;
    double green = (cos(orien) + 1.0)/2.0;
    double blue = 1.0 - red - green;

    red *= mag*4.0;
    green *= mag*4.0;
    blue *= mag*4.0;

    // Make sure the pixel values range from 0 to 255
    red = min(255.0, max(0.0, red));
    green = min(255.0, max(0.0, green));
    blue = min(255.0, max(0.0, blue));

    image->setPixel(c, r, qRgb( (int) (red), (int) (green), (int) (blue)));

    ************************************************************************/
}


void MainWindow::BilinearInterpolation(QImage *image, double x, double y, double rgb[3])
{
    // Add your code here.  Return the RGB values for the pixel at location (x,y) in double rgb[3].
}

// Here is some sample code for rotating an image.  I assume orien is in degrees.

void MainWindow::RotateImage(QImage *image, double orien)
{
    int r, c;
    QRgb pixel;
    QImage buffer;
    int w = image->width();
    int h = image->height();
    double radians = -2.0*3.141*orien/360.0;

    buffer = image->copy();

    pixel = qRgb(0, 0, 0);
    image->fill(pixel);

    for(r=0;r<h;r++)
    {
        for(c=0;c<w;c++)
        {
            double rgb[3];
            double x0, y0;
            double x1, y1;

            // Rotate around the center of the image.
            x0 = (double) (c - w/2);
            y0 = (double) (r - h/2);

            // Rotate using rotation matrix
            x1 = x0*cos(radians) - y0*sin(radians);
            y1 = x0*sin(radians) + y0*cos(radians);

            x1 += (double) (w/2);
            y1 += (double) (h/2);

            BilinearInterpolation(&buffer, x1, y1, rgb);

            image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
        }
    }

}

void MainWindow::FindPeaksImage(QImage *image, double thres)
{
    // Add your code here.
}


void MainWindow::MedianImage(QImage *image, int radius)
{
    // Add your code here
}

void MainWindow::HoughImage(QImage *image)
{
    // Add your code here
}

void MainWindow::CrazyImage(QImage *image)
{
    // Add your code here
}

void MainWindow::RandomSeedImage(QImage *image, int num_clusters)
{
     // Add your code here
}

void MainWindow::PixelSeedImage(QImage *image, int num_clusters)
{
    // Add your code here
}

void MainWindow::HistogramSeedImage(QImage *image, int num_clusters)
{
    // Add your code here
}


// void Convolve(QImage *image, int radius, double* kernel, int height, int width) {
//     int r, c, rd, cd, i;

//     QRgb pixel;

//     // get width and height of image
//     int w = image->width();
//     int h = image->height();

//     int size = 2 * radius + 1;

//     QImage buffer;

//     // make copy of original image with 
//     buffer = image->copy(-radius, -radius, w + 2 * radius, h + 2 * radius);

//     // convolve image with kernel
//     for (r = 0; r < h; r++) {
//         for (c = 0; c < w; c++) {
//             double rgb[3];
//             rgb[0] = 0.0;
//             rgb[1] = 0.0;
//             rgb[2] = 0.0;

//             // for each pixel convolve
//             for (rd = -radius; rd < -radius + height; rd++) {
//                 for (cd = -radius; cd < -radius + width; cd++) {
//                     pixel = buffer.pixel(c + cd + radius, r + rd + radius);
//                     double weight = kernel[width * (rd + radius) + cd + radius]

//                     rgb[0] += weight * (double) qRed(pixel);
//                     rgb[1] += weight * (double) qGreen(pixel);
//                     rgb[2] += weight * (double) qBlue(pixel);
//                 }
//             }

//             image->setPixel(c, r, qRgb((int) floor(rgb[0] + 0.5), (int) floor(rgb[1] + 0.5), (int) floor(rgb[2] + 0.5)));
//         }
//     }
// }