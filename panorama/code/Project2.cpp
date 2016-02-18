#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>
#include "Matrix.h"

/*******************************************************************************
    The following are helper routines with code already written.
    The routines you'll need to write for the assignment are below.
*******************************************************************************/

/*******************************************************************************
Draw detected Harris corners
    interestPts - interest points
    numInterestsPts - number of interest points
    imageDisplay - image used for drawing

    Draws a red cross on top of detected corners
*******************************************************************************/
void MainWindow::DrawInterestPoints(CIntPt *interestPts, int numInterestsPts, QImage &imageDisplay)
{
   int i;
   int r, c, rd, cd;
   int w = imageDisplay.width();
   int h = imageDisplay.height();

   for(i=0;i<numInterestsPts;i++)
   {
       c = (int) interestPts[i].m_X;
       r = (int) interestPts[i].m_Y;

       for(rd=-2;rd<=2;rd++)
           if(r+rd >= 0 && r+rd < h && c >= 0 && c < w)
               imageDisplay.setPixel(c, r + rd, qRgb(255, 0, 0));

       for(cd=-2;cd<=2;cd++)
           if(r >= 0 && r < h && c + cd >= 0 && c + cd < w)
               imageDisplay.setPixel(c + cd, r, qRgb(255, 0, 0));
   }
}

/*******************************************************************************
Compute interest point descriptors
    image - input image
    interestPts - array of interest points
    numInterestsPts - number of interest points

    If the descriptor cannot be computed, i.e. it's too close to the boundary of
    the image, its descriptor length will be set to 0.

    I've implemented a very simple 8 dimensional descriptor.  Feel free to
    improve upon this.
*******************************************************************************/
void MainWindow::ComputeDescriptors(QImage image, CIntPt *interestPts, int numInterestsPts)
{
    int r, c, cd, rd, i, j;
    int w = image.width();
    int h = image.height();
    double *buffer = new double [w*h];
    QRgb pixel;

    // Descriptor parameters
    double sigma = 2.0;
    int rad = 4;

    // Computer descriptors from green channel
    for(r=0;r<h;r++)
       for(c=0;c<w;c++)
        {
            pixel = image.pixel(c, r);
            buffer[r*w + c] = (double) qGreen(pixel);
        }

    // Blur
    SeparableGaussianBlurImage(buffer, w, h, sigma);

    // Compute the desciptor from the difference between the point sampled at its center
    // and eight points sampled around it.
    for(i=0;i<numInterestsPts;i++)
    {
        int c = (int) interestPts[i].m_X;
        int r = (int) interestPts[i].m_Y;

        if(c >= rad && c < w - rad && r >= rad && r < h - rad)
        {
            double centerValue = buffer[(r)*w + c];
            int j = 0;

            for(rd=-1;rd<=1;rd++)
                for(cd=-1;cd<=1;cd++)
                    if(rd != 0 || cd != 0)
                {
                    interestPts[i].m_Desc[j] = buffer[(r + rd*rad)*w + c + cd*rad] - centerValue;
                    j++;
                }

            interestPts[i].m_DescSize = DESC_SIZE;
        }
        else
        {
            interestPts[i].m_DescSize = 0;
        }
    }

    delete [] buffer;
}

/*******************************************************************************
Draw matches between images
    matches - matching points
    numMatches - number of matching points
    image1Display - image to draw matches
    image2Display - image to draw matches

    Draws a green line between matches
*******************************************************************************/
void MainWindow::DrawMatches(CMatches *matches, int numMatches, QImage &image1Display, QImage &image2Display)
{
    int i;
    // Show matches on image
    QPainter painter;
    painter.begin(&image1Display);
    QColor green(0, 250, 0);
    QColor red(250, 0, 0);

    for(i=0;i<numMatches;i++)
    {
        painter.setPen(green);
        painter.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter.setPen(red);
        painter.drawEllipse((int) matches[i].m_X1-1, (int) matches[i].m_Y1-1, 3, 3);
    }

    QPainter painter2;
    painter2.begin(&image2Display);
    painter2.setPen(green);

    for(i=0;i<numMatches;i++)
    {
        painter2.setPen(green);
        painter2.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter2.setPen(red);
        painter2.drawEllipse((int) matches[i].m_X2-1, (int) matches[i].m_Y2-1, 3, 3);
    }

}


/*******************************************************************************
Given a set of matches computes the "best fitting" homography
    matches - matching points
    numMatches - number of matching points
    h - returned homography
    isForward - direction of the projection (true = image1 -> image2, false = image2 -> image1)
*******************************************************************************/
bool MainWindow::ComputeHomography(CMatches *matches, int numMatches, double h[3][3], bool isForward)
{
    int error;
    int nEq=numMatches*2;

    dmat M=newdmat(0,nEq,0,7,&error);
    dmat a=newdmat(0,7,0,0,&error);
    dmat b=newdmat(0,nEq,0,0,&error);

    double x0, y0, x1, y1;

    for (int i=0;i<nEq/2;i++)
    {
        if(isForward == false)
        {
            x0 = matches[i].m_X1;
            y0 = matches[i].m_Y1;
            x1 = matches[i].m_X2;
            y1 = matches[i].m_Y2;
        }
        else
        {
            x0 = matches[i].m_X2;
            y0 = matches[i].m_Y2;
            x1 = matches[i].m_X1;
            y1 = matches[i].m_Y1;
        }


        //Eq 1 for corrpoint
        M.el[i*2][0]=x1;
        M.el[i*2][1]=y1;
        M.el[i*2][2]=1;
        M.el[i*2][3]=0;
        M.el[i*2][4]=0;
        M.el[i*2][5]=0;
        M.el[i*2][6]=(x1*x0*-1);
        M.el[i*2][7]=(y1*x0*-1);

        b.el[i*2][0]=x0;
        //Eq 2 for corrpoint
        M.el[i*2+1][0]=0;
        M.el[i*2+1][1]=0;
        M.el[i*2+1][2]=0;
        M.el[i*2+1][3]=x1;
        M.el[i*2+1][4]=y1;
        M.el[i*2+1][5]=1;
        M.el[i*2+1][6]=(x1*y0*-1);
        M.el[i*2+1][7]=(y1*y0*-1);

        b.el[i*2+1][0]=y0;

    }
    int ret=solve_system (M,a,b);
    if (ret!=0)
    {
        freemat(M);
        freemat(a);
        freemat(b);

        return false;
    }
    else
    {
        h[0][0]= a.el[0][0];
        h[0][1]= a.el[1][0];
        h[0][2]= a.el[2][0];

        h[1][0]= a.el[3][0];
        h[1][1]= a.el[4][0];
        h[1][2]= a.el[5][0];

        h[2][0]= a.el[6][0];
        h[2][1]= a.el[7][0];
        h[2][2]= 1;
    }

    freemat(M);
    freemat(a);
    freemat(b);

    return true;
}


/*******************************************************************************
*******************************************************************************
*******************************************************************************

    The routines you need to implement are below

*******************************************************************************
*******************************************************************************
*******************************************************************************/


/*******************************************************************************
Blur a single channel floating point image with a Gaussian.
    image - input and output image
    w - image width
    h - image height
    sigma - standard deviation of Gaussian

    This code should be very similar to the code you wrote for assignment 1.
*******************************************************************************/
void MainWindow::SeparableGaussianBlurImage(double *image, int w, int h, double sigma)
{
    if(sigma == 0) return;

    int r, c, rd, cd, i;
    double pixel;

    // This is the size of the kernel
    int isig = (int)sigma;
    int size = 2 * isig + 1;


    // Create a buffer image so we're not reading and writing to the same image during filtering.
    double *buffer = new double [w * h];
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            buffer[r * w + c] = image[r * w + c]; 
        }
    }

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
            double result = 0.0;

            // Convolve the kernel at each pixel
            for(rd = -isig; rd <= isig; rd++) {
                if(r + rd >= 0 && r + rd < h) {
            
                    pixel = buffer[(r + rd) * w + c];

                    // Get the value of the kernel
                    double weight = kernel[rd + isig];

                    result += weight * pixel;
                }

            }

            // Store convolved pixel in the image to be returned.
            image[r * w + c] = result;
        }
    }

    // make a new copy from updated image
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            buffer[r * w + c] = image[r * w + c]; 
        }
    }

    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            double result = 0.0;

            // Convolve the kernel at each pixel
            for(cd = -isig; cd <= isig; cd++) {
                if(c + cd >= 0 && c + cd < w) {
                    pixel = buffer[r * w + c + cd];

                    // Get the value of the kernel
                    double weight = kernel[cd + isig];

                    result += weight * pixel;
                }
            }

            // Store convolved pixel in the image to be returned.
            image[r * w + c] = result;
        }
    }

    // Clean up.
    delete [] kernel;
    delete [] buffer;

}


/*******************************************************************************
Detect Harris corners.
    image - input image
    sigma - standard deviation of Gaussian used to blur corner detector
    thres - Threshold for detecting corners
    interestPts - returned interest points
    numInterestsPts - number of interest points returned
    imageDisplay - image returned to display (for debugging)
*******************************************************************************/
void MainWindow::HarrisCornerDetector(QImage image, double sigma, double thres, CIntPt **interestPts, int &numInterestsPts, QImage &imageDisplay)
{
    int r, c;
    int w = image.width();
    int h = image.height();
    double *buffer = new double [w*h];
    QRgb pixel;

    numInterestsPts = 0;

    // Set up matrix
    double *ix = new double[w * h];
    double *iy = new double[w * h];
    double *ixy = new double[w * h];
    double *harris = new double[w * h];

    // Compute the corner response using green band only
    // Initialize ixix, iyiy, ixiy
    for(r = 0; r < h; r++) {
        for(c = 0; c < w; c++) {
            pixel = image.pixel(c, r);

            buffer[r * w + c] = (double) qGreen(pixel);

            ix[r * w + c] = 0.0;
            iy[r * w + c] = 0.0;
            ixy[r * w + c] = 0.0;
        }
    }

    

    // Compute the ixix, iyiy, ixiy
    for(r = 1; r < h - 1; r++) {
        for (c = 1; c < w - 1; c++) {
            double left = buffer[r * w + c - 1];
            double right = buffer[r * w + c + 1];
            double bottom = buffer[(r + 1) * w + c];
            double top = buffer[(r - 1) * w + c];

            ix[r * w + c] = pow(right - left, 2);
            iy[r * w + c] = pow(bottom - top, 2);
            ixy[r * w + c] = (right - left) * (bottom - top);
        }
    }

    // Smooth
    SeparableGaussianBlurImage(ix, w, h, sigma);
    SeparableGaussianBlurImage(iy, w, h, sigma);
    SeparableGaussianBlurImage(ixy, w, h, sigma);

    // Compute Harris = det(H) / trace(H)
    for(r = 1; r < h - 1; r++) {
        for(c = 1; c < w - 1; c++) {
            double xx = ix[r * w + c];
            double yy = iy[r * w + c];
            double xy = ixy[r * w + c];

            double det = xx * yy - xy * xy;
            double trace = xx + yy;

            double val = 0.0;
            if(trace > 0.0)
                val = det / trace;

            harris[r * w + c] = val;

            // Scale for displaying
            val = max(50.0, min(255.0, val));

            imageDisplay.setPixel(c, r, qRgb(val, val, val));
        }
    }

    // Compute interest points and display
    for(r = 1; r < h - 1; r++) {
        for(c = 1; c < w - 1; c++) {

            int x, y;
            bool max = true;
            // The corner response > threshold
            // Might be interesting to look at
            if(harris[r * w + c] > thres) {
                for(x = -1; x <= 1; x++) {
                    for(y = -1; y <= 1; y++) {
                        // Not local miximum of surrounding pixels
                        if(harris[r * w + c] < harris[(r + x) * w + c + y]) {
                            max = false;
                        }

                    }
                }
                // Local maximum, corner
                if(max) {
                    numInterestsPts++;
                    imageDisplay.setPixel(c, r, qRgb((int)255, (int)0, (int)0));
                }
            }
        }
    }


    // Store interest points
    // The descriptor of the interest point is stored in m_Desc
    // The length of the descriptor is m_DescSize, if m_DescSize = 0, then it is not valid.
    *interestPts = new CIntPt [numInterestsPts];
    int i = 0;
    for(r = 1; r < h - 1; r++) {
        for(c = 1; c < w - 1; c++) {
            pixel = imageDisplay.pixel(c, r);
            // Interest point value
            if(qRed(pixel) == 255 && qGreen(pixel) == 0 && qBlue(pixel) == 0) {
                (*interestPts)[i].m_X = (double) c;
                (*interestPts)[i].m_Y = (double) r;
                i++;
            }
        }
    }


    // Once you are done finding the interest points, display them on the image
    DrawInterestPoints(*interestPts, numInterestsPts, imageDisplay);

    // Clean up
    delete [] buffer;
    delete [] ix;
    delete [] ixy;
    delete [] iy;
    delete [] harris;
}


/*******************************************************************************
Find matching interest points between images.
    image1 - first input image
    interestPts1 - interest points corresponding to image 1
    numInterestsPts1 - number of interest points in image 1
    image2 - second input image
    interestPts2 - interest points corresponding to image 2
    numInterestsPts2 - number of interest points in image 2
    matches - set of matching points to be returned
    numMatches - number of matching points returned
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::MatchInterestPoints(QImage image1, CIntPt *interestPts1, int numInterestsPts1,
                             QImage image2, CIntPt *interestPts2, int numInterestsPts2,
                             CMatches **matches, int &numMatches, QImage &image1Display, QImage &image2Display)
{

    // Compute the descriptors for each interest point.
    // You can access the descriptor for each interest point using interestPts1[i].m_Desc[j].
    // If interestPts1[i].m_DescSize = 0, it was not able to compute a descriptor for that point
    ComputeDescriptors(image1, interestPts1, numInterestsPts1);
    ComputeDescriptors(image2, interestPts2, numInterestsPts2);

    
    *matches = new CMatches[numInterestsPts1];

    int x, y, z;
    numMatches = 0;

    // Compute best match from image 2 for each interest point in image 1
    // Best matching point has the smallest norm distance
    for(x = 0; x < numInterestsPts1; x++) {
        // Check if descriptor valid
        if(interestPts1[x].m_DescSize > 0) {
            CIntPt pt1 = interestPts1[x];
            int min = 0;
            double mindist = INFINITY;

            for(y = 0; y < numInterestsPts2; y++) {
                // Check if descriptor valid
                if(interestPts2[y].m_DescSize > 0) {
                    CIntPt pt2 = interestPts2[y];
                    double dist = 0.0;
                    for(z = 0; z < DESC_SIZE; z++) {
                        dist += pow(pt1.m_Desc[z] - pt2.m_Desc[z], 2);
                    }
                    if(dist < mindist) {
                        mindist = dist;
                        min = y;
                    }
                }
            }

            // Store matching points
            (*matches)[numMatches].m_X1 = interestPts1[x].m_X;
            (*matches)[numMatches].m_Y1 = interestPts1[x].m_Y;
            (*matches)[numMatches].m_X2 = interestPts2[min].m_X;
            (*matches)[numMatches].m_Y2 = interestPts2[min].m_Y;
            numMatches++;
        }
    }


    // The position of the interest point in iamge 1 is (m_X1, m_Y1)

    // Draw the matches
    DrawMatches(*matches, numMatches, image1Display, image2Display);
}

/*******************************************************************************
Project a point (x1, y1) using the homography transformation h
    (x1, y1) - input point
    (x2, y2) - returned point
    h - input homography used to project point
*******************************************************************************/
void MainWindow::Project(double x1, double y1, double &x2, double &y2, double h[3][3])
{
    double u = h[0][0] * x1 + h[0][1] * y1 + h[0][2];
    double v = h[1][0] * x1 + h[1][1] * y1 + h[1][2];
    double w = h[2][0] * x1 + h[2][1] * y1 + h[2][2];

    x2 = u / w;
    y2 = v / w;
}

/*******************************************************************************
Count the number of inliers given a homography.  This is a helper function for RANSAC.
    h - input homography used to project points (image1 -> image2
    matches - array of matching points
    numMatches - number of matchs in the array
    inlierThreshold - maximum distance between points that are considered to be inliers

    Returns the total number of inliers.
*******************************************************************************/
int MainWindow::ComputeInlierCount(double h[3][3], CMatches *matches, int numMatches, double inlierThreshold)
{
    int i;
    int count = 0;

    for(i = 0; i < numMatches; i++) {
        double x1 = matches[i].m_X1;
        double y1 = matches[i].m_Y1;
        double x2, y2;

        // Project the first point
        Project(x1, y1, x2, y2, h);

        // Check projected point an inlier if distance < threshold
        double distance = pow(x2 - matches[i].m_X2, 2) + pow(y2 - matches[i].m_Y2, 2);

        if(distance < inlierThreshold * inlierThreshold) count++;

    }

    return count;
}


/*******************************************************************************
Compute homography transformation between images using RANSAC.
    matches - set of matching points between images
    numMatches - number of matching points
    numIterations - number of iterations to run RANSAC
    inlierThreshold - maximum distance between points that are considered to be inliers
    hom - returned homography transformation (image1 -> image2)
    homInv - returned inverse homography transformation (image2 -> image1)
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::RANSAC(CMatches *matches, int numMatches, int numIterations, double inlierThreshold,
                        double hom[3][3], double homInv[3][3], QImage &image1Display, QImage &image2Display)
{
    int i, j;
    CMatches randomMatches[4];
    double h[3][3];

    int max = 0;

    for(i = 0; i < numIterations; i++) {

        // Randomly select 4 matching pairs
        
        for(j = 0; j < 4; j++) {
            int num = rand() % numMatches;

            randomMatches[j].m_X1 = matches[num].m_X1;
            randomMatches[j].m_X2 = matches[num].m_X2;
            randomMatches[j].m_Y1 = matches[num].m_Y1;
            randomMatches[j].m_Y2 = matches[num].m_Y2;
        }

        // Compute Homography for the 4 random matches and inliers count
        if(ComputeHomography(randomMatches, 4, h, true)) {
            int count = ComputeInlierCount(h, matches, numMatches, inlierThreshold);

            // Check if the homography is the best
            if(count > max) {
                max = count;
                int x, y;
                for(x = 0; x < 3; x++) {
                    for(y = 0; y < 3; y++) {
                        hom[x][y] = h[x][y];
                    }
                }
            }
        }
    }

    CMatches *inliers = new CMatches[max];
    int numInliers = 0;

    for(i = 0; i < numMatches; i++) {
        double x2, y2;

        Project(matches[i].m_X1, matches[i].m_Y1, x2, y2, hom);

        double distance = pow(x2 - matches[i].m_X2, 2.0) + pow(y2 - matches[i].m_Y2, 2.0);

        if(distance < inlierThreshold * inlierThreshold) {
            inliers[numInliers].m_X1 = matches[i].m_X1;
            inliers[numInliers].m_Y1 = matches[i].m_Y1;
            inliers[numInliers].m_X2 = matches[i].m_X2;
            inliers[numInliers].m_Y2 = matches[i].m_Y2;

            numInliers++;
        }

    }

    // Recompute the best homography and inverse homography
    ComputeHomography(inliers, numInliers, hom, true);
    ComputeHomography(inliers, numInliers, homInv, false);

    // After you're done computing the inliers, display the corresponding matches.
    DrawMatches(inliers, numInliers, image1Display, image2Display);

    // Clean up
    delete [] inliers;

}

/*******************************************************************************
Bilinearly interpolate image (helper function for Stitch)
    image - input image
    (x, y) - location to interpolate
    rgb - returned color values

    You can just copy code from previous assignment.
*******************************************************************************/
bool MainWindow::BilinearInterpolation(QImage *image, double x, double y, double rgb[3])
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
    if(x1 >= 0 &&  x1 < w && x2 >= 0 && x2 < w && y1 >= 0 && y1 < h && y2 >= 0 && y2 < h) {
        QRgb ltop = image->pixel(x1, y1);
        QRgb rtop = image->pixel(x1, y2);
        QRgb lbot = image->pixel(x2, y1);
        QRgb rbot = image->pixel(x2, y2);

        rgb[0] += (x2 - x) * (y2 - y) * (double)qRed(ltop) + (x2 - x) * (y - y1) * (double)qRed(rtop) + 
                (x - x1) * (y2 - y) * (double)qRed(lbot) + (x - x1) * (y - y1) * (double)qRed(rbot);
        rgb[1] += (x2 - x) * (y2 - y) * (double)qGreen(ltop) + (x2 - x) * (y - y1) * (double)qGreen(rtop) + 
                (x - x1) * (y2 - y) * (double)qGreen(lbot) + (x - x1) * (y - y1) * (double)qGreen(rbot);
        rgb[2] += (x2 - x) * (y2 - y) * (double)qBlue(ltop) + (x2 - x) * (y - y1) * (double)qBlue(rtop) + 
                (x - x1) * (y2 - y) * (double)qBlue(lbot) + (x - x1) * (y - y1) * (double)qBlue(rbot);

    }
    else
        return false;

    return true;
}


/*******************************************************************************
Stitch together two images using the homography transformation
    image1 - first input image
    image2 - second input image
    hom - homography transformation (image1 -> image2)
    homInv - inverse homography transformation (image2 -> image1)
    stitchedImage - returned stitched image
*******************************************************************************/
void MainWindow::Stitch(QImage image1, QImage image2, double hom[3][3], double homInv[3][3], QImage &stitchedImage)
{
    int r, c, i;
    QRgb pixel;
    int ws, hs;
    int w1 = image1.width();
    int h1 = image1.height();
    int w2 = image2.width();
    int h2 = image2.height();
    int cmin = 0;
    int cmax = w1;
    int rmin = 0;
    int rmax = h1;
    double corners[8];

    corners[0] = corners[1] = corners[3] = corners[4] = 0.0;
    corners[2] = corners[6] = (double) w2;
    corners[5] = corners[7] = (double) h2;

    // Compute ws and has by projecting four corners of image2 to image1
    for(i = 0; i < 4; i++) {
        double x, y;
        Project(corners[i * 2], corners[i * 2 + 1], x, y, homInv);
        int r = floor(y + 0.5);
        int c = floor(x + 0.5);
        cmin = cmin > c ? c : cmin;
        cmax = cmax < c ? c : cmax;
        rmin = rmin > r ? r : rmin;
        rmax = rmax < r ? r : rmax;
    }

    ws = cmax - cmin;
    hs = rmax - rmin;

    // Initialize the stitchedImage
    stitchedImage = QImage(ws, hs, QImage::Format_RGB32);
    stitchedImage.fill(qRgb(0,0,0));

    // Warp image1 and image2 to stitchedImage.
    // Copy image1 to stitchedImage
    for(r = 0; r < h1; r++) {
        for(c = 0; c < w1; c++) {
            pixel = image1.pixel(c, r);
            stitchedImage.setPixel(c - cmin, r - rmin, pixel);
        }
    }

    // Project each pixel of stitchedImage onto image2
    // bilinearInterpolate the value
    for(r = 0; r < hs; r++) {
        for(c = 0; c < ws; c++) {
            double x, y; 
            double rgb[3];
            pixel = stitchedImage.pixel(c, r);
            Project((double)(c + cmin), (double)(r + rmin), x, y, hom);
            if(BilinearInterpolation(&image2, x, y, rgb)) {
                stitchedImage.setPixel(c, r, qRgb((int) rgb[0], (int) rgb[1], (int) rgb[2]));
            }

        }
    }

}

