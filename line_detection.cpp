/*`
 line_detection.cpp
 polygon_detection
 
 Created by Andrew Ho on 6/24/18.
 Copyright © 2018 Andrew Ho, California Institute of Technology. All rights reserved.
 */

#include "line_detection.hpp"

void restrictEdges(vector<gNode*> filtered_corners)
{
    
    for(auto && tpoint : filtered_corners)
    {
        if(tpoint->getNumNeighbors() > 2)
        {
            float d1 = 0;
            float d2 = 0;
            
            gNode* sig_1;
            gNode* sig_2;
            
            for(auto && neighbor : tpoint->adjacent)
            {
                float test_dist = calcMagnitudeGNode(tpoint, neighbor);
                if(test_dist > d1)
                {
                    d1 = test_dist;
                    sig_1 = neighbor;
                }
                else if (test_dist < d1 && test_dist > d2)
                {
                    d2 = test_dist;
                    sig_2 = neighbor;
                }
            }
            
            for(auto && neighbor : tpoint->adjacent)
            {
                if(neighbor != sig_1 && neighbor != sig_2)
                {
                    tpoint->removeEdge(neighbor);
                    neighbor->removeEdge(tpoint);
                }
            }
        }
    }
}

void filterEdges(vector<gNode*> filtered_corners){
    for(auto && tpoint : filtered_corners)
    {
        Point anchor(tpoint->x, tpoint->y);
        int n1 = 0;
        while(n1 < tpoint->getNumNeighbors()-1)
        {
            gNode* gNode1 = tpoint->adjacent[n1];
            Point p1(tpoint->adjacent[n1]->x, tpoint->adjacent[n1]->y);
            int n2 = 1;
            while(n2 < tpoint->getNumNeighbors())
            {
                gNode* gNode2 = tpoint->adjacent[n2];
                Point p2(tpoint->adjacent[n2]->x, tpoint->adjacent[n2]->y);
                float angle = angleBetween(p1, anchor, p2);
                
                if(angle < 20)
                {
                    if(calcMagnitude(p1, anchor) > calcMagnitude(anchor, p2))
                    {
                        tpoint->removeEdge(gNode2);
                        gNode2->removeEdge(tpoint);
                        
                    }
                    else
                    {
                        tpoint->removeEdge(gNode1);
                        gNode2->removeEdge(tpoint);
                        break;
                    }
                }
                else
                {
                    ++n2;
                }
            }
            ++n1;
        }
    }
    
    for(auto && tpoint : filtered_corners)
    {
        Point anchor(tpoint->x, tpoint->y);
        int n1 = 0;
        while(n1 < tpoint->getNumNeighbors())
        {
            gNode* gNode1 = tpoint->adjacent[n1];
            Point p1(tpoint->adjacent[n1]->x, tpoint->adjacent[n1]->y);
            
            int n2 = 0;
            while(n2 < gNode1->getNumNeighbors())
            {
                gNode* gNode2 = tpoint->adjacent[n2];
                Point p2(gNode1->adjacent[n2]->x, gNode1->adjacent[n2]->y);
                
                float angle = angleBetween(anchor, p1, p2);
                
                if(angle > 160)
                {
                    tpoint->removeEdge(gNode1);
                    gNode1->removeEdge(tpoint);
                    
                    if(!(tpoint->edgeExist(gNode2))){
                        tpoint->addEdge(gNode2);
                        gNode2->addEdge(tpoint);
                    }
                    break;
                }
                else
                {
                    ++n2;
                }
            }
            ++n1;
        }
    }
}



bool testPixels(Mat *frame, int test_x, int test_y, int window_size){
    int low = (-1) * (window_size / 2);
    int high = window_size + low;
    bool pixExist = false;
    for(int i = low; i < high; ++i){
        for(int j = low; j < high; ++j){
            if ((i != 0) && (j != 0)){
                int pix=(int)frame->at<uchar>(test_y+i, test_x+j);
                if(pix != 0){
                    pixExist = true;
                    break;
                }
            }
        }
    }
    
    return pixExist;
}

/**
    Takes an image and performs morphological operations (i.e opening/closing)
    repeatedly
 
 @param input_image : The image to perform operations on
 @param operation : The type of morphological operation to be performed
 @param iter : The number of times to perform the specified operation
 
 */
void morphOp(Mat *input_image, int operation, int iter)
{
    // Perform morphological closing
    if(operation == 0)
    {
        for(int i = 0; i < iter; ++i)
        {
            erode(*input_image, *input_image, Mat());
            dilate(*input_image, *input_image, Mat());
        }
    }
    // Perform morphological opening
    else
    {
        for(int i = 0; i < iter; ++i)
        {
            dilate(*input_image, *input_image, Mat());
            erode(*input_image, *input_image, Mat());
        }
    }
}

/**
 Takes a set of points and adds edges between them by checking for non-zero pixels
 
 @param frame : The scaffold used to check for pixels
 @param filtered_corners : The type of morphological operation to be performed
 
 */
void findEdges(Mat frame, vector<gNode*> filtered_corners)
{
    for(int i = 0; i < filtered_corners.size() - 1; ++i)
    {
        gNode* p1 = filtered_corners[i];
        for(int j = i + 1; j < filtered_corners.size(); ++j)
        {
            gNode* p2 = filtered_corners[j];
            
            // If the edge does not exist, then check pixels between
            
            // Step parameters
            float step = 0.05;
            int num_steps = 1/step + 1;
            int WPix_threshold = static_cast<int>(0.7 * num_steps);
            
            // Initialize WPix_count to zero
            int WPix_count = 0;
            int empty_stretch = 0;
            // Test pixels
            for(int param_t = 1; param_t < num_steps; param_t++)
            {
                float test_x = p1->x + (param_t*step)*(p2->x - p1->x);
                float test_y = p1->y + (param_t*step)*(p2->y - p1->y);
                
                bool pixExist = testPixels(&frame, test_x, test_y, 5);
                
                
                if(pixExist)
                {
                    step = 0.05;
                    ++WPix_count;
                    empty_stretch = 0;
                }
                else{
                    step = 0.02;
                    empty_stretch++;
                    if(empty_stretch > 3){
                        break;
                    }
                }
            }
            if(WPix_count > WPix_threshold)
            {
                if(calcMagnitude(Point(p1->x, p1->y), Point(p2->x, p2->y)) > 20)
                {
                    p1->addEdge(p2);
                    p2->addEdge(p1);
                }

            }
        }
    }
}

/**
 Finds the intersection between three sets of detected points
 
 @param ST_points : Points detected by Shi-Tomasi corner detection
 @param Hough_points : Points detected by Hough Line Transform
 @param Hough_tree : K-D Tree using HLT points, used for knn search
 @param LSD_points : Points detected by Line Segment Detector
 @param LSD_tree : K-D Tree using LSD points, used for knn search
 
 @return filtered_corners : Intersection between ST, HLT, and LSD points
 */
vector<gNode> filterPoints(vector<Point2f> ST_points,
                           cv::Mat_<float> Hough_points,
                           cv::flann::Index *Hough_tree,
                           cv::Mat_<float> LSD_points,
                           cv::flann::Index *LSD_tree)
{
    
    vector<gNode> filtered_corners;
    
    for(auto && point : ST_points)
    {
        unsigned int num_neighbours = 1;
        
        Mat query = (cv::Mat_<float>(1, 2) << point.x, point.y);
        
        Mat hough_indices, hough_dists;
        Mat LSD_indices, LSD_dists;
        
        Hough_tree->knnSearch(query, hough_indices, hough_dists, num_neighbours, cv::flann::SearchParams(32));
        
        LSD_tree->knnSearch(query, LSD_indices, LSD_dists, num_neighbours, cv::flann::SearchParams(32));
        
        float hough_dist = hough_dists.at<float>(0,0);
        float LSD_dist =  LSD_dists.at<float>(0,0);
        
        if(hough_dist <= 75 && LSD_dist <= 75)
        {
            Point2f hough_pt(Hough_points.at<float>(hough_indices.at<int>(0,0), 0), Hough_points.at<float>(hough_indices.at<int>(0,0), 1));
            Point2f LSD_pt(LSD_points.at<float>(LSD_indices.at<int>(0,0), 0), LSD_points.at<float>(LSD_indices.at<int>(0,0), 1));
            Point2f median = calcCenter(point, hough_pt, LSD_pt);
            
            gNode medianGNode(median);
            
            filtered_corners.push_back(medianGNode);
        }
    }
    
    return filtered_corners;
}

/**
 Takes in a set of points and converts them to the correct dimensions for
 K-D Tree construction
 
 @param input_lines : Vector of detected points
 
 @return flannPoints : Converted points
 */
cv::Mat_<float> convertPoints(vector<Vec4f> input_lines)
{
    cv::Mat_<float> flannPoints(0,2);

    for(auto && line : input_lines)
    {
        cv::Mat p1 = (cv::Mat_<float>(1, 2) << line[0], line[1]); //First endpoint
        cv::Mat p2 = (cv::Mat_<float>(1, 2) << line[2], line[3]); //Second endpoint
        flannPoints.push_back(p1);
        flannPoints.push_back(p2);
    }
    
    return flannPoints;
}

/**
 Takes set of points and builds a K-D Tree
 
 @param flannPoints : Points to be added to K-D Tree
 
 @return KDtree : Constructed K-D Tree
 */
cv::flann::Index createKDTree(cv::Mat_<float> flannPoints)
{
    
    cv::flann::Index KDtree(flannPoints, cv::flann::KDTreeIndexParams(1));
 
    return KDtree;
}

/**
 Crops image to the bounding box
 
 @param src : Original source image to be cropped
 @param ROI : Points detected by Hough Line Transform
 @param backgroundPresent : Environment Variable to indicate which pathway
 
 @return cropped_WGE : Cropped image
 */
Mat cropImg(Mat src, Rect ROI, bool backgroundPresent)
{
    Mat cropped_WGE;
    if(backgroundPresent)
    {
        Mat prewitt = prewittOperator(src, 1.0);
        cvtColor(prewitt, prewitt, COLOR_BGR2GRAY);
        Mat quantile_99 = WGE(prewitt, 99.5);
        cropped_WGE = crop(quantile_99, ROI);
        
    }
    else
    {
        cropped_WGE = crop(src, ROI);
        cvtColor(cropped_WGE, cropped_WGE, COLOR_BGR2GRAY);
    }
    
    return cropped_WGE;
}

/**
 Helper function to crop image for a specified ROI
 
 @param input_img : Source image to be cropped
 @param ROI : Rectangle of Interest
 
 @return cropped_WGE : Cropped image
 */
Mat crop(Mat input_img, Rect ROI)
{
    return input_img(ROI);
}

/**
 Applies heavy gradient elimintation to determine rectangle of interest containing the satellite
 
 @param input : Source image
 @param backgroundPresent : Environment variable, used to determine pipeline
 
 @return ROI : Rectangle dimensions and coordinates containing the tightest bounding box around the satellite
 */
Rect findROI(Mat input, bool backgroundPresent)
{
    // Find tightest ROI
    Mat ROI_temp;
    
    Mat testing_img;
    
    // Prepare testing image
    if(backgroundPresent)
    {
        Mat prewitt = prewittOperator(input, 1.0);
        
        // Convert image to grayscale
        cvtColor(input, input, COLOR_BGR2GRAY);
        
        // Blur image for ROI search
        cv::Mat prewitt_blur;
        GaussianBlur(prewitt, prewitt_blur, Size(5,5), 0, 0);
        
        // Convert image to grayscale
        cvtColor(prewitt_blur, prewitt_blur, COLOR_BGR2GRAY);
        
        testing_img = WGE(prewitt_blur, 99.95);
    }
    else
    {
        cvtColor(input, input, COLOR_BGR2GRAY);
        testing_img = input;
    }
    
    ROI_temp = adaptive_thresh(testing_img, 11, 12, 1);
    
    int leftBound = input.cols;
    int rightBound = 0;
    bool topFlag = false;
    int topBound = input.rows;
    int botBound = 0;
    
    for(int height_iter = 0; height_iter < input.rows; ++height_iter)
    {
        for(int width_iter = 0; width_iter < input.cols; ++width_iter)
        {
            int pix=(int)ROI_temp.at<uchar>(height_iter, width_iter);
            
            if(pix != 0)
            {
                if(!topFlag)
                {
                    topFlag = true;
                    topBound = height_iter;
                }
                else
                {
                    if(height_iter > botBound)
                    {
                        botBound = height_iter;
                    }
                }
                
                if(width_iter < leftBound)
                {
                    leftBound = width_iter;
                }
                else if (width_iter > rightBound)
                {
                    rightBound = width_iter;
                }
            }
        }
    }
    
    Rect ROI(leftBound, topBound, rightBound - leftBound, botBound - topBound);
    
    // Resize ROI
    int diff = abs(ROI.width - ROI.height);
    int buffer = max(ROI.width, ROI.height) * 0.3;
    
    if(ROI.height > ROI.width)
    {
        ROI.x = ROI.x - (diff/2);
        ROI.width = ROI.width + diff;
    }
    else
    {
        ROI.y = ROI.y - (diff/2);
        ROI.height = ROI.height + diff;
    }
    
    ROI.x = max(0, ROI.x - buffer);
    ROI.y = max(0, ROI.y - buffer);
    
    ROI.width = min(ROI.width + 2*buffer, input.cols-ROI.x);
    ROI.height = min(ROI.height + 2*buffer, input.rows-ROI.y);
    
    ROI.width = min(ROI.width, ROI.height);
    ROI.height = min(ROI.width, ROI.height);
    
    return ROI;
}

/**
 Calculates histogram of input image and removes pixels below the specified threshold
 
 @param input_img : Source image
 @param threshold : Pixel value threshold
 
 @return removedGradient : Image after weak gradient pixels are removed
 */
Mat WGE(Mat input_img, float threshold)
{
    // calculate histogram for every pixel value (i.e [0 - 255])
    cv::Mat hist;
    int histSize = 256;
    float range[] = { 0, 256 } ;
    const float* histRange = { range };
    bool uniform = true; bool accumulate = false;
    cv::calcHist( &input_img, 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, uniform, accumulate );
    
    // total pixels in image
    float totalPixels = input_img.cols * input_img.rows;
    
    // calculate percentage of every histogram bin (i.e: pixel value [0 - 255])
    // the 'bins' variable holds pairs of (int pixelV alue, float percentage)
    std::vector<std::pair<int, float>> bins;
    float percentage;
    for(int i = 0; i < 256; ++i)
    {
        percentage = (hist.at<float>(i,0)*100.0)/totalPixels;
        bins.push_back(std::make_pair(i, percentage));
    }
    
    // sort the bins according to percentage
    sort(bins.begin(), bins.end());
    
    // compute percentile for a pixel value
    
    Mat removedGradient = input_img.clone();
    
    for (int i = 0; i < input_img.rows; i++)
    {
        for (int j = 0; j < input_img.cols; j++)
        {
            int pixel = input_img.at<uchar>(i,j);
            
            float sum = 0;
            for (auto b : bins){
                if(b.first != pixel)
                {
                    sum += b.second;
                }
                else
                {
                    sum += b.second/2;
                    break;
                }
            }
            
            if(sum < threshold)
            {
                removedGradient.at<uchar>(i,j) = 0;
            }
        }
    }
    
    return removedGradient;
}

/**
 Helper function to find the directional gradient at a pixel
 
 @param input_image : Source image to be cropped
 @param x : X coordinate of pixel
 @param y : Y coordinate of pixel
 
 @return standardized direction vector with float values
 */

Vec2f calcPixGradient(Mat input_image, int x, int y)
{
    // Calculate vertical gradient
    float v_grad;

    int pixel1 = input_image.at<Vec3b>(y-1, x-1)[0] * -3;
    int pixel2 = input_image.at<Vec3b>(y, x-1)[0] * 0;
    int pixel3 = input_image.at<Vec3b>(y+1, x-1)[0] * 3;
    
    int pixel4 = input_image.at<Vec3b>(y-1, x)[0] * -10;
    int pixel5 = input_image.at<Vec3b>(y, x)[0] * 0;
    int pixel6 = input_image.at<Vec3b>(y+1, x)[0] * 10;
    
    int pixel7 = input_image.at<Vec3b>(y-1, x+1)[0] * -3;
    int pixel8 = input_image.at<Vec3b>(y, x+1)[0] * 0;
    int pixel9 = input_image.at<Vec3b>(y+1, x+1)[0] * 3;
    
    int sum = abs(pixel1 + pixel2 + pixel3 + pixel4 + pixel5 + pixel6 + pixel7 + pixel8 + pixel9);
    v_grad = min(sum, 255);
    
    
    // Calculate horizontal gradient
    float h_grad;
    
    pixel1 = input_image.at<Vec3b>(y-1, x-1)[0] * -3;
    pixel2 = input_image.at<Vec3b>(y, x-1)[0] * -10;
    pixel3 = input_image.at<Vec3b>(y+1, x-1)[0] * -3;
    
    pixel4 = input_image.at<Vec3b>(y-1, x)[0] * 0;
    pixel5 = input_image.at<Vec3b>(y, x)[0] * 0;
    pixel6 = input_image.at<Vec3b>(y+1, x)[0] * 0;
    
    pixel7 = input_image.at<Vec3b>(y-1, x+1)[0] * 3;
    pixel8 = input_image.at<Vec3b>(y, x+1)[0] * 10;
    pixel9 = input_image.at<Vec3b>(y+1, x+1)[0] * 3;
    
    sum = abs(pixel1 + pixel2 + pixel3 + pixel4 + pixel5 + pixel6 + pixel7 + pixel8 + pixel9);
    h_grad = min(sum, 255);
    
    // Standardize vectors
    float norm = sqrt(h_grad*h_grad+v_grad*v_grad) ;
    
    return Vec2f(h_grad/norm, v_grad/norm);
}

/**
 Helper function to apply convolution filter with specified coefficient
 
 @param input_image : Image to apply filter to
 @param coeff : Convlutional coefficient
 
 @return finalPrewitt : Final product of convolutional filters
 */
Mat prewittOperator(Mat input_image, float coeff)
{
    Mat verticalPrewitt = input_image.clone();
    for (int i = 0; i < input_image.rows; i++)
    {
        for (int j = 0; j < input_image.cols; j++)
        {
            if(i == 0 || i == input_image.rows || j == 0 || j == input_image.cols)
            {
                verticalPrewitt.at<Vec3b>(i,j)[0] = 0;
                verticalPrewitt.at<Vec3b>(i,j)[1] = 0;
                verticalPrewitt.at<Vec3b>(i,j)[2] = 0;
            }
            else
            {
                int pixel1 = input_image.at<Vec3b>(i-1,j-1)[0] * -1;
                int pixel2 = input_image.at<Vec3b>(i,j-1)[0] * 0;
                int pixel3 = input_image.at<Vec3b>(i+1,j-1)[0] * 1;
                
                int pixel4 = input_image.at<Vec3b>(i-1,j)[0] * -1 * coeff;
                int pixel5 = input_image.at<Vec3b>(i,j)[0] * 0;
                int pixel6 = input_image.at<Vec3b>(i+1,j)[0] * 1 * coeff;
                
                int pixel7 = input_image.at<Vec3b>(i-1,j+1)[0] * -1;
                int pixel8 = input_image.at<Vec3b>(i,j+1)[0] * 0;
                int pixel9 = input_image.at<Vec3b>(i+1,j+1)[0] * 1;
                
                int sum = abs(pixel1 + pixel2 + pixel3 + pixel4 + pixel5 + pixel6 + pixel7 + pixel8 + pixel9);
                sum = min(sum, 255);
                
                verticalPrewitt.at<Vec3b>(i,j)[0] = sum;
                verticalPrewitt.at<Vec3b>(i,j)[1] = sum;
                verticalPrewitt.at<Vec3b>(i,j)[2] = sum;
            }
        }
    }

    Mat horizontalPrewitt = input_image.clone();
    for (int i = 0; i < input_image.rows; i++)
    {
        for (int j = 0; j < input_image.cols; j++)
        {
            if(i == 0 || i == input_image.rows || j == 0 || j == input_image.cols)
            {
                horizontalPrewitt.at<Vec3b>(i, j)[0] = 0;
                horizontalPrewitt.at<Vec3b>(i, j)[1] = 0;
                horizontalPrewitt.at<Vec3b>(i, j)[2] = 0;
            }
            else
            {
                int pixel1 = input_image.at<Vec3b>(i-1,j-1)[0] * -1;
                int pixel2 = input_image.at<Vec3b>(i,j-1)[0] * -1 * coeff;
                int pixel3 = input_image.at<Vec3b>(i+1,j-1)[0] * -1;
            
                int pixel4 = input_image.at<Vec3b>(i-1,j)[0] * 0;
                int pixel5 = input_image.at<Vec3b>(i,j)[0] * 0;
                int pixel6 = input_image.at<Vec3b>(i+1,j)[0] * 0;
            
                int pixel7 = input_image.at<Vec3b>(i-1,j+1)[0] * 1;
                int pixel8 = input_image.at<Vec3b>(i,j+1)[0] * 1 * coeff;
                int pixel9 = input_image.at<Vec3b>(i+1,j+1)[0] * 1;
            
                int sum = abs(pixel1 + pixel2 + pixel3 + pixel4 + pixel5 + pixel6 + pixel7 + pixel8 + pixel9);
                sum = min(sum, 255);
            
                horizontalPrewitt.at<Vec3b>(i,j)[0] = sum;
                horizontalPrewitt.at<Vec3b>(i,j)[1] = sum;
                horizontalPrewitt.at<Vec3b>(i,j)[2] = sum;
            }
        }
    }
    
    Mat finalPrewitt = verticalPrewitt + horizontalPrewitt;
    
    return finalPrewitt;
}

/**
 Calculates centroid of three specified points
 
 @param p1 : First point
 @param p2 : Second point
 @param p3 : Third point
 
 @return retPoint : Point with the coordinates of the centroid
 */
Point2f calcCenter(Point2f p1, Point2f p2, Point2f p3)
{
    Point2f retPoint;
    
    retPoint.x = (p1.x + p2.x + p3.x)/3;
    retPoint.y = (p1.y + p2.y + p3.y)/3;
    
    return retPoint;
}

/**
 Takes input image and performs one of five thresholding operations for a specified threshold
 
 @param src_gray : Source image to be thresholded
 @param thresh_Val : Threshold value
 
 @return threshold_output : Result image from thresholding
 */

Mat global_thresh(Mat src_gray, int thresh_val, int input_operation)
{
        int threshold_value = thresh_val;
        int const max_binary_value = 255;
    
    
        /* 0: Binary
         1: Binary Inverted
         2: Threshold Truncated
         3: Threshold to Zero
         4: Threshold to Zero Inverted
         */
        int operation = input_operation;
    
        Mat threshold_output;
    
        threshold(src_gray,
                  threshold_output,
                  threshold_value,
                  max_binary_value,
                  operation);
    
        return threshold_output;
}

/**
 Calculates Otsu's value using histogram analysis
 
 @param src : Source image
 
 @return max_val : Otsu's value
 */
double getOtsuVal( const cv::Mat& _src )
{
    cv::Size size = _src.size();
    if ( _src.isContinuous() )
    {
        size.width *= size.height;
        size.height = 1;
    }
    const int N = 256;
    int i, j, h[N] = {0};
    for ( i = 0; i < size.height; i++ )
    {
        const uchar* src = _src.data + _src.step*i;
        for ( j = 0; j <= size.width - 4; j += 4 )
        {
            int v0 = src[j], v1 = src[j+1];
            h[v0]++; h[v1]++;
            v0 = src[j+2]; v1 = src[j+3];
            h[v0]++; h[v1]++;
        }
        for ( ; j < size.width; j++ )
        {
            h[src[j]]++;
        }
    }

    double mu = 0, scale = 1./(size.width*size.height);
    for ( i = 0; i < N; i++ )
    {
        mu += i*h[i];
    }
    
    mu *= scale;
    double mu1 = 0, q1 = 0;
    double max_sigma = 0, max_val = 0;

    for ( i = 0; i < N; i++ )
    {
        double p_i, q2, mu2, sigma;

        p_i = h[i]*scale;
        mu1 *= q1;
        q1 += p_i;
        q2 = 1. - q1;

        if ( std::min(q1,q2) < FLT_EPSILON || std::max(q1,q2) > 1. - FLT_EPSILON )
        {
            continue;
        }
        
        mu1 = (mu1 + i*p_i)/q1;
        mu2 = (mu - q1*mu1)/q2;
        sigma = q1*q2*(mu1 - mu2)*(mu1 - mu2);
        if ( sigma > max_sigma )
        {
            max_sigma = sigma;
            max_val = i;
        }
    }

    return max_val;
}

/**
 Thinning algorithm to find 1 pixel wide skeleton
 
 @param im : image to thin
 @param iter : number of iterations to perform thinning
 */

void thinningIteration(cv::Mat& im, int iter)
{
    cv::Mat marker = cv::Mat::zeros(im.size(), CV_8UC1);
    
    for (int i = 1; i < im.rows-1; i++)
    {
        for (int j = 1; j < im.cols-1; j++)
        {
            uchar p2 = im.at<uchar>(i-1, j);
            uchar p3 = im.at<uchar>(i-1, j+1);
            uchar p4 = im.at<uchar>(i, j+1);
            uchar p5 = im.at<uchar>(i+1, j+1);
            uchar p6 = im.at<uchar>(i+1, j);
            uchar p7 = im.at<uchar>(i+1, j-1);
            uchar p8 = im.at<uchar>(i, j-1);
            uchar p9 = im.at<uchar>(i-1, j-1);
            
            int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
            (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
            (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
            (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);
            int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
            int m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8);
            int m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8);
            
            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
            {
                marker.at<uchar>(i,j) = 1;
            }
        }
    }
    
    im &= ~marker;
}

/**
 Helper function for thinning algorithm
 
 @param src : Source image to be cropped
 @param ROI : Rectangle of Interest
 
 */
void thinning(cv::Mat& im)
{
    im /= 255;
    
    cv::Mat prev = cv::Mat::zeros(im.size(), CV_8UC1);
    cv::Mat diff;
    
    do
    {
        thinningIteration(im, 0);
        thinningIteration(im, 1);
        cv::absdiff(im, prev, diff);
        im.copyTo(prev);
    }
    while (cv::countNonZero(diff) > 0);
    
    im *= 255;
}

/**
 Performs adaptive thresholding on a specified image
 
 @param src_gray : Input image to threshold
 @param input_blocksize : Blocksize, dictates window size for adaptive thresholding
 @param input_cval : Constant weight to readjust window values
 @param input_operation : Type of thresholding to perform
 
 @return threshold_output : Image after adaptive thresholding is performed
 */

Mat adaptive_thresh(Mat src_gray, int input_blocksize, int input_cval, int input_operation)
{
    int blocksize = input_blocksize;
    int c_val = input_cval;

    /* 0: Binary
     1: Binary Inverted
     2: Threshold Truncated
     3: Threshold to Zero
     4: Threshold to Zero Inverted
     */
    int operation = input_operation;

    Mat threshold_output;

    adaptiveThreshold(src_gray,
                      threshold_output,
                      255,
                      ADAPTIVE_THRESH_MEAN_C,
                      operation,
                      blocksize,
                      c_val);

    return threshold_output;
}

/**
 Helper function to calculate the distance between two points
 
 @param p1 : First point
 @param p2 : Second point
 
 @return Distance between input points
 */
float calcMagnitude(Point p1, Point p2)
{
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

float calcMagnitudeGNode(gNode* p1, gNode* p2)
{
    return sqrt(pow(p2->x - p1->x, 2) + pow(p2->y - p1->y, 2));
}

/**
 Helper function to calculate dot product between two vectors
 
 @param v1 : First vector
 @param v2 : Second vector
 
 @return Dot product between the input vectors
 */
float dotproduct(vector<float> v1, vector<float> v2)
{
    return (v1[0]*v2[0]) + (v1[1]*v2[1]);
}

/**
 Helper function to calculate angle between three points
 Points are specified like so:
 
 p1
  \
   \
    p2 ------p3
 
 @param p1 : First point
 @param p2 : Second point
 @param p3 : Third point
 
 @return angle_deg : Angle between the three points in degrees
 */
float angleBetween(Point p1, Point p2, Point p3)
{
    vector<float> v1;
    v1.push_back(p1.x - p2.x);
    v1.push_back(p1.y - p2.y);
    
    vector<float> v2;
    v2.push_back(p3.x - p2.x);
    v2.push_back(p3.y - p2.y);
    
    float dot_product = dotproduct(v1, v2);
    
    float mag1 = calcMagnitude(p1, p2);
    float mag2 = calcMagnitude(p3, p2);
    
    float angle_rad = acos(dot_product/(mag1 * mag2));
    
    float angle_deg = angle_rad * (180/3.14);
    
    return angle_deg;
}

/**
 Performs Shi-Tomasi Corner detection, detects a specified amount of corners
 
 @param maxCorners : Number of corners to detect
 @param input : Input image to find corners in
 
 @return corners : Vector of corner coordinates
 */
vector<Point2f> goodFeaturesToTrack_Callback(int maxCorners, Mat input )
{
    /// Parameters for Shi-Tomasi algorithm
    maxCorners = MAX(maxCorners, 1);
    vector<Point2f> corners;
    double qualityLevel = 0.01;
    double minDistance = 10;
    int blockSize = 21, gradientSize = 21;
    bool useHarrisDetector = true;
    double k = 0.04;
    
    /// Apply corner detection
    goodFeaturesToTrack(input,
                        corners,
                        maxCorners,
                        qualityLevel,
                        minDistance,
                        Mat(),
                        blockSize,
                        gradientSize,
                        useHarrisDetector,
                        k );
    
    return corners;
}
