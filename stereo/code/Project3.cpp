#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>

/*******************************************************************************
    The following are helper routines with code already written.
    The routines you'll need to write for the assignment are below.
*******************************************************************************/

/*******************************************************************************
   K-means segment the image
        image - Input image
        gridSize - Initial size of the segments
        numIterations - Number of iterations to run k-means
        spatialSigma - Spatial sigma for measuring distance
        colorSigma - Color sigma for measuring distance
        matchCost - The match cost for each pixel at each disparity
        numDisparities - Number of disparity levels
        segmentImage - Image showing segmentations

*******************************************************************************/
void MainWindow::Segment(QImage image, int gridSize, int numIterations, double spatialSigma, double colorSigma,
                         double *matchCost, int numDisparities, QImage *segmentImage)
{
    int w = image.width();
    int h = image.height();
    int iter;
    int numSegments = 0;

    // Stores the segment assignment for each pixel
    int *segment = new int [w*h];

    // Compute an initial segmentation
    GridSegmentation(segment, numSegments, gridSize, w, h);

    // allocate memory for storing the segments mean position and color
    double (*meanSpatial)[2] = new double [numSegments][2];
    double (*meanColor)[3] = new double [numSegments][3];

    // Iteratively update the segmentation
    for(iter=1;iter<numIterations;iter++)
    {
        // Compute new means
        ComputeSegmentMeans(image, segment, numSegments, meanSpatial, meanColor);
        // Compute new pixel assignment to pixels
        AssignPixelsToSegments(image, segment, numSegments, meanSpatial, meanColor, spatialSigma, colorSigma);
    }

    // Update means again for display
    ComputeSegmentMeans(image, segment, numSegments, meanSpatial, meanColor);
    // Display the segmentation
    DrawSegments(segmentImage, segment, meanColor);

    // Update the match cost based on the segmentation
    SegmentAverageMatchCost(segment, numSegments, w, h, numDisparities, matchCost);

    delete [] meanSpatial;
    delete [] meanColor;
    delete [] segment;
}

/*******************************************************************************
   Compute initial segmentation of the image using a grid
        segment - Segment assigned to each pixel
        numSegments - Number of segments
        gridSize - Size of the grid-based segments
        w - Image width
        h - Image height

*******************************************************************************/
void MainWindow::GridSegmentation(int *segment, int &numSegments, int gridSize, int w, int h)
{
    int r, c;
    int step = w/gridSize;

    if(step*gridSize < w)
        step += 1;

    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
             int rs = r/gridSize;
             int cs = c/gridSize;

             segment[r*w + c] = rs*step + cs;

             numSegments = rs*step + cs + 1;

        }

}

/*******************************************************************************
   Draw the image segmentation
        segmentImage - Image to display the segmentation
        segment - Segment assigned to each pixel
        meanColor - The mean color of the segments

*******************************************************************************/
void MainWindow::DrawSegments(QImage *segmentImage, int *segment, double (*meanColor)[3])
{
    int w = segmentImage->width();
    int h = segmentImage->height();
    int r, c;

    for(r=0;r<h-1;r++)
        for(c=0;c<w-1;c++)
        {
            int segIdx = segment[r*w + c];
            if(segIdx != segment[r*w + c + 1] ||
               segIdx != segment[(r+1)*w + c])
            {
                segmentImage->setPixel(c, r, qRgb(255, 255, 255));
            }
            else
            {
                segmentImage->setPixel(c, r, qRgb((int) meanColor[segIdx][0],
                                                  (int) meanColor[segIdx][1], (int) meanColor[segIdx][2]));
            }
        }
}

/*******************************************************************************
   Display the computed disparities
        disparities - The disparity for each pixel
        disparityScale - The amount to scale the disparity for display
        minDisparity - Minimum disparity
        disparityImage - Image to display the disparity
        errorImage - Image to display the error
        GTImage - The ground truth disparities
        m_DisparityError - The average error

*******************************************************************************/
void MainWindow::DisplayDisparities(double *disparities, int disparityScale, int minDisparity,
                        QImage *disparityImage, QImage *errorImage, QImage GTImage, double *disparityError)
{
    int w = disparityImage->width();
    int h = disparityImage->height();
    int r, c;
    int gtw = GTImage.width();
    bool useGT = false;
    double pixelCt = 0.0;
    *disparityError = 0.0;
    double maxError = 1.0*(double) disparityScale;

    if(gtw == w)
        useGT = true;

    QRgb pixel;

    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
            double disparity = disparities[r*w + c];
            disparity *= (double) disparityScale;
            disparity -= minDisparity*disparityScale;

            disparityImage->setPixel(c, r, qRgb((int) disparity, (int) disparity, (int) disparity));

            if(useGT)
            {
                pixel = GTImage.pixel(c, r);

                if(qGreen(pixel) > 0)
                {
                    double dist = fabs(disparity - (double) qGreen(pixel));
                    if(dist > maxError)
                        (*disparityError)++;
                    pixelCt++;

                    if(dist > maxError)
                        errorImage->setPixel(c, r, qRgb(255,255,255));
                    else
                        errorImage->setPixel(c, r, qRgb(0,0,0));
                }


            }
        }

    if(useGT)
        *disparityError /= pixelCt;
}

/*******************************************************************************
   Render warped views between the images
        image - Image to be warped
        disparities - The disparities for each pixel
        disparityScale - The amount to warp the image, usually between 0 and 1
        renderImage - The final rendered image

*******************************************************************************/
void MainWindow::Render(QImage image, double *disparities, double disparityScale, QImage *renderImage)
{
    int r, c;
    int w = image.width();
    int h = image.height();
    double *projDisparity = new double [w*h];
    double *projDisparityCt = new double [w*h];
    QRgb pixel0;
    QRgb pixel1;

    memset(projDisparity, 0, w*h*sizeof(double));
    memset(projDisparityCt, 0, w*h*sizeof(double));

    // First forward project the disparity values
    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
            double disparity =  -disparities[r*w + c]*disparityScale;
            double x = (double) c + disparity;
            int cp = (int) x;
            double del = x - (double) cp;

            if(cp >= 0 && cp < w-1)
            {
                // Make sure we get the depth ordering correct.
                if(projDisparityCt[r*w + cp] == 0)
                {
                    projDisparity[r*w + cp] = (1.0 - del)*disparity;
                    projDisparityCt[r*w + cp] = (1.0 - del);
                }
                else
                {
                    // Make sure the depth ordering is correct
                    if(fabs(disparity) > fabs(2.0 + projDisparity[r*w + cp]/projDisparityCt[r*w + cp]))
                    {
                        projDisparity[r*w + cp] = (1.0 - del)*disparity;
                        projDisparityCt[r*w + cp] = (1.0 - del);
                    }
                    else
                    {
                        projDisparity[r*w + cp] += (1.0 - del)*disparity;
                        projDisparityCt[r*w + cp] += (1.0 - del);
                    }
                }

                if(projDisparityCt[r*w + cp + 1] == 0)
                {
                    projDisparity[r*w + cp + 1] = (del)*disparity;
                    projDisparityCt[r*w + cp + 1] = (del);
                }
                else
                {
                    // Make sure the depth ordering is correct
                    if(fabs(disparity) > fabs(2.0 + projDisparity[r*w + cp + 1]/projDisparityCt[r*w + cp + 1]))
                    {
                        projDisparity[r*w + cp + 1] = (del)*disparity;
                        projDisparityCt[r*w + cp + 1] = (del);
                    }
                    else
                    {
                        projDisparity[r*w + cp + 1] += (del)*disparity;
                        projDisparityCt[r*w + cp + 1] += (del);
                    }
                }
            }
        }

    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
        {
            if(projDisparityCt[r*w + c] > 0.0)
            {
                projDisparity[r*w + c] /= projDisparityCt[r*w + c];
            }
        }

    // Fill in small holes after the forward projection
    FillHoles(projDisparity, projDisparityCt, w, h);

    renderImage->fill(qRgb(0,0,0));

    // Backward project to find the color values for each pixel
    for(r=0;r<h;r++)
        for(c=0;c<w;c++)
            if(projDisparityCt[r*w + c] > 0.0)
        {
            double disparity =  projDisparity[r*w + c];
            double x = (double) c - disparity;
            int cp = (int) x;
            double del = x - (double) cp;

            if(cp >= 0 && cp < w-1)
            {
                pixel0 = image.pixel(cp, r);
                pixel1 = image.pixel(cp+1, r);

                int red = (int) ((1.0 - del)*(double)qRed(pixel0) + del*(double)qRed(pixel1));
                int green = (int) ((1.0 - del)*(double)qGreen(pixel0) + del*(double)qGreen(pixel1));
                int blue = (int) ((1.0 - del)*(double)qBlue(pixel0) + del*(double)qBlue(pixel1));

                // Uncomment if you want to see the disparities
            //    red = (int) disparity*4.0;
            //    green = (int) disparity*4.0;
            //    blue = (int) disparity*4.0;

                renderImage->setPixel(c, r, qRgb(red, green, blue));
            }
        }


    delete [] projDisparity;
    delete [] projDisparityCt;
}

/*******************************************************************************
   Fill holes in the projected disparities (Render helper function)
        projDisparity - Projected disparity
        projDisparityCt - The weight of each projected disparity.  A value of 0 means the pixel doesn't have a disparity
        w, h - The width and height of the image

*******************************************************************************/
void MainWindow::FillHoles(double *projDisparity, double *projDisparityCt, int w, int h)
{
    int r, c, cd, rd;
    double *bufferCt = new double [w*h];

    memcpy(bufferCt, projDisparityCt, w*h*sizeof(double));

    for(r=1;r<h-1;r++)
        for(c=1;c<w-1;c++)
            if(bufferCt[r*w + c] == 0)
        {
            double avgDisparity = 0.0;
            double avgCt = 0.0;

            for(rd=-1;rd<=1;rd++)
                for(cd=-1;cd<=1;cd++)
                {
                    int idx = (r + rd)*w + c + cd;
                   avgDisparity += projDisparity[idx]*bufferCt[idx];
                   avgCt += bufferCt[idx];
                }

            if(avgCt > 0.0)
            {
                projDisparity[r*w + c] = avgDisparity/avgCt;
                projDisparityCt[r*w + c] = avgCt;

            }
        }

    delete [] bufferCt;
}


/*******************************************************************************
*******************************************************************************
*******************************************************************************

    The routines you need to implement are below

*******************************************************************************
*******************************************************************************
*******************************************************************************/

/*******************************************************************************
Compute match cost using Squared Distance
    image1 - Input image 1
    image2 - Input image 2
    minDisparity - Minimum disparity between image 1 and image 2
    maxDisparity - Minimum disparity between image 1 and image 2
    offset - offset of window to compute the SSD score inside the window of size [2*offset+1]
    matchCost - The match cost (squared distance) between pixels

    To access the match cost at pixel (c, r) at disparity d use
    matchCost[d*w*h + r*w + c]
*******************************************************************************/
void MainWindow::SSD(QImage image1, QImage image2, int minDisparity, int maxDisparity, int offset, double *matchCost)
{
	// SSD is performed by sum(image1's pixel - image2's pixel)^2
    int w = image1.width();
	int h = image1.height();

    int r, c, d, i, j;
    double rd, gd, bd;
    QRgb p1;
    QRgb p2;

    for(d = minDisparity; d < maxDisparity; d++) {
        for(r = offset; r < h - offset; r++) {
            for(c = offset; c < w - offset; c++) {
            	matchCost[(d - minDisparity) * w * h + r * w + c] = 0.0;
                if(c - d >= offset && c - d < w - offset) {
                    
                    rd = gd = bd = 0.0; // aggregated sum over window

                    for(i = -offset; i <= offset; i++) {
                        for(j = -offset; j <= offset; j++) {
                            p1 = image1.pixel(c + i, r + j);
                            p2 = image2.pixel(c - d + i, r + j);

							rd += pow((double)qRed(p1) - (double)qRed(p2), 2.0);
		    				gd += pow((double)qGreen(p1) - (double)qGreen(p2), 2.0);
		    				bd += pow((double)qBlue(p1) - (double)qBlue(p2), 2.0);

                        }
                    }
                    rd /= pow(offset * 2 + 1, 2.0);
                    gd /= pow(offset * 2 + 1, 2.0);
                    bd /= pow(offset * 2 + 1, 2.0);
                    matchCost[(d - minDisparity) * w * h + r * w + c] = sqrt(rd + gd + bd);
                } 
            }
        }
    }


}

/*******************************************************************************
Compute match cost using Absolute Distance
    image1 - Input image 1
    image2 - Input image 2
    minDisparity - Minimum disparity between image 1 and image 2
    maxDisparity - Minimum disparity between image 1 and image 2
    offset - offset of window to compute the SAD score inside the window of size [2*offset+1]
    matchCost - The match cost (absolute distance) between pixels

    To access the match cost at pixel (c, r) at disparity d use
    matchCost[d*w*h + r*w + c]
*******************************************************************************/
void MainWindow::SAD(QImage image1, QImage image2, int minDisparity, int maxDisparity, int offset, double *matchCost)
{
	// SAD is performed by sum(abs(image1's pixel - image2's pixel))
    int w = image1.width();
	int h = image1.height();

    int r, c, d, i, j;
    double rd, gd, bd;
    QRgb p1;
    QRgb p2;

    for(d = minDisparity; d < maxDisparity; d++) {

        for(r = offset; r < h - offset; r++) {
            for(c = offset; c < w - offset; c++) {
            	matchCost[(d - minDisparity) * w * h + r * w + c] = 0.0;
                if(c - d >= offset && c - d < w - offset) {
                    
                    rd = gd = bd = 0.0; // aggregated sum over window

                    for(i = -offset; i <= offset; i++) {
                        for(j = -offset; j <= offset; j++) {
                            p1 = image1.pixel(c + i, r + j);
                            p2 = image2.pixel(c - d + i, r + j);

                         
							rd += abs((double)qRed(p1) - (double)qRed(p2));
		    				gd += abs((double)qGreen(p1) - (double)qGreen(p2));
		    				bd += abs((double)qBlue(p1) - (double)qBlue(p2));

                        }
                    }
                    rd /= pow(offset * 2 + 1, 2.0);
                    gd /= pow(offset * 2 + 1, 2.0);
                    bd /= pow(offset * 2 + 1, 2.0);
                    matchCost[(d - minDisparity) * w * h + r * w + c] = rd + gd + bd;
                } 
            }
        }
    }
}

/*******************************************************************************
Compute match cost using Normalized Cross Correlation
    image1 - Input image 1
    image2 - Input image 2
    minDisparity - Minimum disparity between image 1 and image 2
    maxDisparity - Minimum disparity between image 1 and image 2
    offset - offset of window to compute the NCC score inside the window of size [2*offset+1]
    matchCost - The match cost (1 - NCC) between pixels

    To access the match cost at pixel (c, r) at disparity d use
    matchCost[d*w*h + r*w + c]
*******************************************************************************/
void MainWindow::NCC(QImage image1, QImage image2, int minDisparity, int maxDisparity, int offset, double *matchCost)
{
	// NCC is performed by sum(image1 * image2) / sqrt(sum(image1)^2 * sum(image)^2)
	int w = image1.width();
	int h = image1.height();

    int r, c, d, i, j;
    double i1, i2, i12;
    QRgb p1;
    QRgb p2;

    for(d = minDisparity; d < maxDisparity; d++) {

        for(r = offset; r < h - offset; r++) {
            for(c = offset; c < w - offset; c++) {
            	matchCost[(d - minDisparity) * w * h + r * w + c] = 0.0;
                if(c - d >= offset && c - d < w - offset) {
                    
                    i1 = i2 = i12 = 0.0;

                    for(i = -offset; i <= offset; i++) {
                        for(j = -offset; j <= offset; j++) {
                            p1 = image1.pixel(c + i, r + j);
                            p2 = image2.pixel(c - d + i, r + j);

                            double r1 = (double) qRed(p1);
                            double g1 = (double) qGreen(p1);
                            double b1 = (double) qBlue(p1);

                            double r2 = (double) qRed(p2);
                            double g2 = (double) qGreen(p2);
                            double b2 = (double) qBlue(p2);

                            i1 += r1 * r1 + g1 * g1 + b1 * b1;

                            i12 += r1 * r2 + g1 * g2 + b1 * b2;

                            i2 += r2 * r2 + g2 * g2 + b2 * b2;
                        }
                    }

                    matchCost[(d - minDisparity) * w * h + r * w + c] = 1.0 - i12 / (sqrt(i1 * i2));

                } 
            }
        }
    }

}

/*******************************************************************************
Gaussian blur the match score.
    matchCost - The match cost between pixels
    w, h - The width and height of the image
    numDisparities - The number of disparity levels
    sigma - The standard deviation of the blur kernel

    I would recommend using SeparableGaussianBlurImage as a helper function.
*******************************************************************************/
void MainWindow::GaussianBlurMatchScore(double *matchCost, int w, int h, int numDisparities, double sigma)
{
    int d;

    for(d = 0; d < numDisparities; d++) {
        SeparableGaussianBlurImage(&(matchCost[d * w * h]), w, h, sigma);
    }
}

/*******************************************************************************
Blur a floating piont image using Gaussian kernel (helper function for GaussianBlurMatchScore.)
    image - Floating point image
    w, h - The width and height of the image
    sigma - The standard deviation of the blur kernel

    You may just cut and paste code from previous assignment
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
Bilaterally blur the match score using the colorImage to compute kernel weights
    matchCost - The match cost between pixels
    w, h - The width and height of the image
    numDisparities - The number of disparity levels
    sigmaS, sigmaI - The standard deviation of the blur kernel for spatial and intensity
    colorImage - The color image
*******************************************************************************/
void MainWindow::BilateralBlurMatchScore(double *matchCost, int w, int h, int numDisparities,
                                         double sigmaS, double sigmaI, QImage colorImage)
{
     double  *buffer;

    QRgb pixel;

    // Size an dradius of the kernel
    int radius = (int)(ceil(3 * sigmaS));
    int size = ((2 * radius) + 1);


    buffer = new double [w*h*numDisparities];
    memcpy(buffer, matchCost, w*h*numDisparities*sizeof(double));
    memset(matchCost, 0, w*h*numDisparities*sizeof(double));

    if(sigmaS <= 0)
        return;

    double *kernel = new double [size];

    for(int i=0;i<size;i++)
    {
        double dist = (double) (i - radius);

        kernel[i] = exp(-(dist*dist)/(2.0*sigmaS*sigmaS));
    }

    for(int i =radius;i<h-radius;i++)
    {
        for(int j=radius;j <w-radius;j++)
        {
            double sum = 0.000001;

            pixel = colorImage.pixel(j, i);
            double value1 = (double) qGreen(pixel);

            for(int rd=-radius;rd<=radius;rd++)
                for(int cd=-radius;cd<=radius;cd++)
                {
                     pixel = colorImage.pixel(j+cd, i+rd);
                     double weight = kernel[rd + radius]*kernel[cd + radius];

                     double value2 = qGreen(pixel);

                     weight *= exp(-((value1 - value2)*(value1 - value2))/(2.0*sigmaI*sigmaI));

                     for(int k=0;k<numDisparities;k++)
                        matchCost[k*h*w + i*w + j] += weight*(double) buffer[k*w*h + (i+rd)*w + j + cd];

                     sum += weight;
                }

            for(int k=0;k<numDisparities;k++)
                matchCost[k*h*w + i*w + j] /= sum;
        }
    }


    delete [] buffer;
    delete [] kernel;
}

/*******************************************************************************
Compute the mean color and position for each segment (helper function for Segment.)
    image - Color image
    segment - Image segmentation
    numSegments - Number of segments
    meanSpatial - Mean position of segments
    meanColor - Mean color of segments
*******************************************************************************/
void MainWindow::ComputeSegmentMeans(QImage image, int *segment, int numSegments, double (*meanSpatial)[2], double (*meanColor)[3])
{
	int w = image.width();
	int h = image.width();
    double *mean = new double[numSegments];
    QRgb pixel;
    
}

/*******************************************************************************
Assign each pixel to the closest segment using position and color
    image - Color image
    segment - Image segmentation
    numSegments - Number of segments
    meanSpatial - Mean position of segments
    meanColor - Mean color of segments
    spatialSigma - Assumed standard deviation of the spatial distribution of pixels in segment
    colorSigma - Assumed standard deviation of the color distribution of pixels in segment
*******************************************************************************/
void MainWindow::AssignPixelsToSegments(QImage image, int *segment, int numSegments, double (*meanSpatial)[2], double (*meanColor)[3],
                            double spatialSigma, double colorSigma)
{
    // Add your code here
}

/*******************************************************************************
Update the match cost based ont eh segmentation.  That is, average the match cost
for each pixel in a segment.
    segment - Image segmentation
    numSegments - Number of segments
    width, height - Width and height of image
    numDisparities - Number of disparities
    matchCost - The match cost between pixels
*******************************************************************************/
void MainWindow::SegmentAverageMatchCost(int *segment, int numSegments,
                                         int w, int h, int numDisparities, double *matchCost)
{
    // Add your code here
}

/*******************************************************************************
For each pixel find the disparity with minimum match cost
    matchCost - The match cost between pixels
    disparities - The disparity for each pixel (use disparity[r*w + c])
    width, height - Width and height of image
    minDisparity - The minimum disparity
    numDisparities - Number of disparities
*******************************************************************************/
void MainWindow::FindBestDisparity(double *matchCost, double *disparities, int w, int h, int minDisparity, int numDisparities)
{
    int r, c, d, i;
    double mincost;

    // Find mincost for in each disparity of each pixel
   
    for(r = 0; r < h; r++) {
    	for(c = 0; c < w; c++) {
    		mincost = INFINITY;
    		i = 0;
    		for(d = 0; d < numDisparities; d++) {
    			if(mincost > matchCost[d * w * h + r * w + c]) {
    				mincost = matchCost[d * w * h + r * w + c];
    				i = d + minDisparity;
    			}
    		}
    		disparities[r * w + c] = (double)i; // store the index
    	}
	}
}

/*******************************************************************************
Create your own "magic" stereo algorithm
    image1 - Input image 1
    image2 - Input image 2
    minDisparity - Minimum disparity between image 1 and image 2
    maxDisparity - Minimum disparity between image 1 and image 2
    param1 - The first parameter to your algorithm
    param2 - The second paramater to your algorithm
    matchCost - The match cost (squared distance) between pixels
*******************************************************************************/
void MainWindow::MagicStereo(QImage image1, QImage image2, int minDisparity, int maxDisparity, double param1, double param2, double *matchCost)
{
    // Add your code here

}
