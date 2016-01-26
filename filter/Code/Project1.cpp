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
    // [-1 0 1] is the filter
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
    // This time is all neighbors - center, [1 -2 1] is the 1d filter
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
    // Graytones then convolve with kernel
    int r, c, rd, cd;

    QRgb pixel;

    // Graytones first
    BlackWhiteImage(image);

    QImage buffer;
    int w = image->width();
    int h = image->height();

    double gx[9] = {1, 0, -1, 2, 0, -2, 1, 0, -1};
    double gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    buffer = image->copy(-1, -1, w + 2, h + 2);

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double magX = 0.0;
            double magY = 0.0;

            for(rd = -1; rd <= 1; rd++) {
                for(cd = -1; cd <= 1; cd++) {
                    pixel = buffer.pixel(c + cd + 1, r + rd + 1);

                    // Get the value of the kernel
                    double weightX = gx[(rd + 1) * 3 + cd + 1];
                    double weightY = gy[(rd + 1) * 3 + cd + 1];

                    magX += weightX * (double) qRed(pixel);
                    magY += weightY * (double) qRed(pixel);
                }
            }

            magX /= 8.0;
            magY /= 8.0;

            double mag = sqrt(magX * magX + magY * magY); // magnitude
            double orien = atan2(magY, magX); // orientation

            double red = (sin(orien) + 1.0) / 2.0;
            double green = (cos(orien) + 1.0) / 2.0;
            double blue = 1.0 - red - green;

            red *= mag * 4.0;
            green *= mag * 4.0;
            blue *= mag * 4.0;

            // Make sure the pixel values range from 0 to 255
            red = min(255.0, max(0.0, red));
            green = min(255.0, max(0.0, green));
            blue = min(255.0, max(0.0, blue));

            // Store convolved pixel in the image to be returned.
            image->setPixel(c, r, qRgb((int)red, (int)green, (int)blue));
        }
    }

}


void MainWindow::BilinearInterpolation(QImage *image, double x, double y, double rgb[3])
{
    // Return the RGB values for the pixel at location (x,y) in double rgb[3].

    rgb[0] = rgb[1] = rgb[2] = 0.0;

    int w = image->width();
    int h = image->height();


    int x1 = static_cast<int>(floor(x));
    int y1 = static_cast<int>(floor(y));
    int x2 = static_cast<int>(ceil(x+0.00001));
    int y2 = static_cast<int>(ceil(y+0.00001));

    // Check if x1, y1, x2, y2 in boundary of image

    QRgb ltop = ((0 <= x1 && x1 < w && 0 <= y1 && y1 < h) ?
                    image->pixel(x1, y1) : qRgb(0, 0, 0));
    QRgb rtop = ((0 <= x1 && x1 < w && 0 <= y2 && y2 < h) ?
                        image->pixel(x1, y2) : qRgb(0, 0, 0));
    QRgb lbot = ((0 <= x2 && x2 < w && 0 <= y1 && y1 < h) ?
                        image->pixel(x2, y1) : qRgb(0, 0, 0));
    QRgb rbot = ((0 <= x2 && x2 < w && 0 <= y2 && y2 < h) ?
                        image->pixel(x2, y2) : qRgb(0, 0, 0));

    double denom = 1/((x2 - x1) * (y2 - y1));
    rgb[0] += (x2 - x) * (y2 - y) * (double)qRed(ltop) + (x2 - x) * (y - y1) * (double)qRed(rtop) + 
            (x - x1) * (y2 - y) * (double)qRed(lbot) + (x - x1) * (y - y1) * (double)qRed(rbot);
    rgb[1] += (x2 - x) * (y2 - y) * (double)qGreen(ltop) + (x2 - x) * (y - y1) * (double)qGreen(rtop) + 
            (x - x1) * (y2 - y) * (double)qGreen(lbot) + (x - x1) * (y - y1) * (double)qGreen(rbot);
    rgb[2] += (x2 - x) * (y2 - y) * (double)qBlue(ltop) + (x2 - x) * (y - y1) * (double)qBlue(rtop) + 
            (x - x1) * (y2 - y) * (double)qBlue(lbot) + (x - x1) * (y - y1) * (double)qBlue(rbot);

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
    // First compute edge magnitude and orientation for the image with sobel operator

    int r, c, rd, cd;

    QRgb pixel;

    // Graytones first
    BlackWhiteImage(image);

    QImage buffer;

    // store pixel magnitude
    QImage magnitude = image->copy();

    // store pixel orientation
    QImage orientation = image->copy();

    int w = image->width();
    int h = image->height();

    double gx[9] = {1, 0, -1, 2, 0, -2, 1, 0, -1};
    double gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

    buffer = image->copy(-1, -1, w + 2, h + 2);

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double magX = 0.0;
            double magY = 0.0;

            for(rd = -1; rd <= 1; rd++) {
                for(cd = -1; cd <= 1; cd++) {
                    pixel = buffer.pixel(c + cd + 1, r + rd + 1);

                    // Get the value of the kernel
                    double weightX = gx[(rd + 1) * 3 + cd + 1];
                    double weightY = gy[(rd + 1) * 3 + cd + 1];

                    magX += weightX * (double) qRed(pixel);
                    magY += weightY * (double) qRed(pixel);
                }
            }

            magX /= 8.0;
            magY /= 8.0;

            double mag = sqrt(magX * magX + magY * magY); // magnitude
            magnitude.setPixel(c, r, qRgb(0, (int)mag, 0));

            double orien = atan2(magY, magX); // orientation

            double red = sin(orien);
            double green = cos(orien);

            orientation.setPixel(c, r, qRgb((int)red, (int)green, 0));

        }
    }

    pixel = qRgb(0, 0, 0);
    // zero out the image
    image->fill(pixel);


    // Compare edge magnitude
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            pixel = magnitude.pixel(c, r);
            double mag0 = (double)qGreen(pixel);

            if(mag0 > thres) {
                pixel = orientation.pixel(c, r);
                double dsin = (double)qRed(pixel);
                double dcos = (double)qGreen(pixel);

                // pixel perpendicular
                double x1, y1;

                double e0[3];

                x1 = dcos + c;
                y1 = dsin + r;

                BilinearInterpolation(&magnitude, x1, y1, e0);

                // pixel perpendicular
                double x2, y2;

                double e1[3];

                x2 = -dcos + c;
                y2 = -dsin + r;

                BilinearInterpolation(&magnitude, x2, y2, e1);

                if(mag0 >= e0[1] && mag0 >= e1[1]) {
                    image->setPixel(c, r, qRgb(255, 255, 255));
                }

            }

        }
    }

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
    int r, c, i, j, k;
    QRgb pixel;
    int iteration = 100;

    int w = image->width();
    int h = image->height();

    // for assign cluster
    int* cluster = new int[h * w];

    // store k means
    int* mean = new int[num_clusters * 3];

    // random select k means
    for(i = 0; i < num_clusters; i++) {
        for(j = 0; j < 3; j++) {
            mean[i * 3 + j] = rand() % 256;
        }
    }

    for(i = 0; i < iteration; i++) {
        for(r = 0; r < h; r++) {
            for(c = 0; c < w; c++) {

                pixel = image->pixel(c, r);
                int min = INFINITY;
                int select = 0;
                for(k = 0; k < num_clusters; k++) {
                    int distance = abs(qRed(pixel) - mean[k * 3]) + abs(qGreen(pixel) - mean[k * 3 + 1]) + abs(qBlue(pixel) - mean[k * 3 + 2]);
                    if (distance < min) {
                        min = distance;
                        select = k;
                    }
                }

                cluster[r * w + c] = select;
            }
        }

        // keep track of pixel count of each cluster
        int* count = new int[num_clusters];
        for(k = 0; k < num_clusters; k++) {
            count[k] = 0;
            mean[k * 3] = mean[k * 3 + 1] = mean[k * 3 + 2] = 0;
            for(r = 0; r < h; r++) {
                for(c = 0; c < w; c++) {
                    pixel = image->pixel(c, r);
                    
                    // accumulate cluster sum
                    if(cluster[r * w + c] == k) {
                        mean[k * 3] += qRed(pixel);
                        mean[k * 3 + 1] += qGreen(pixel);
                        mean[k * 3 + 2] += qBlue(pixel);
                        count[k] += 1;
                    }

                }
            }

            // Average
            if(count[k] != 0) {
                mean[k * 3] /= count[k];
                mean[k * 3 + 1] /= count[k];
                mean[k * 3 + 2] /= count[k];
            }
        }

        delete [] count;
    }

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            image->setPixel(c, r, qRgb(mean[cluster[r * w + c] * 3], 
                mean[cluster[r * w + c] * 3 + 1], mean[cluster[r * w + c] * 3 + 2]));
        }
    }

    delete [] cluster;
    delete [] mean;

}

void MainWindow::PixelSeedImage(QImage *image, int num_clusters)
{
    int r, c, i, j, k;
    QRgb pixel;
    int iteration = 100;

    int w = image->width();
    int h = image->height();

    // for assign cluster
    int* cluster = new int[h * w];

    // store k means
    int* mean = new int[num_clusters * 3];

    // random sample k pixels
    for(i = 0; i < num_clusters; i++) {
        int x = rand() % w;
        int y = rand() % h;
        pixel = image->pixel(x, y);
        mean[i * 3] = qRed(pixel);
        mean[i * 3 + 1] = qGreen(pixel);
        mean[i * 3 + 2] = qBlue(pixel);
    }

    for(i = 0; i < iteration; i++) {
        for(r = 0; r < h; r++) {
            for(c = 0; c < w; c++) {

                pixel = image->pixel(c, r);
                int min = INFINITY;
                int select = 0;
                for(k = 0; k < num_clusters; k++) {
                    int distance = abs(qRed(pixel) - mean[k * 3]) + abs(qGreen(pixel) - mean[k * 3 + 1]) + abs(qBlue(pixel) - mean[k * 3 + 2]);
                    if (distance < min) {
                        min = distance;
                        select = k;
                    }
                }

                cluster[r * w + c] = select;
            }
        }

        // keep track of pixel count of each cluster
        int* count = new int[num_clusters];
        for(k = 0; k < num_clusters; k++) {
            count[k] = 0;
            mean[k * 3] = mean[k * 3 + 1] = mean[k * 3 + 2] = 0;
            for(r = 0; r < h; r++) {
                for(c = 0; c < w; c++) {
                    pixel = image->pixel(c, r);
                    
                    // accumulate cluster sum
                    if(cluster[r * w + c] == k) {
                        mean[k * 3] += qRed(pixel);
                        mean[k * 3 + 1] += qGreen(pixel);
                        mean[k * 3 + 2] += qBlue(pixel);
                        count[k] += 1;
                    }

                }
            }

            // Average
            if(count[k] != 0) {
                mean[k * 3] /= count[k];
                mean[k * 3 + 1] /= count[k];
                mean[k * 3 + 2] /= count[k];
            }
        }

        delete [] count;
    }

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            image->setPixel(c, r, qRgb(mean[cluster[r * w + c] * 3], 
                mean[cluster[r * w + c] * 3 + 1], mean[cluster[r * w + c] * 3 + 2]));
        }
    }

    delete [] cluster;
    delete [] mean;
}

void MainWindow::HistogramSeedImage(QImage *image, int num_clusters)
{
    // Add your code here
}