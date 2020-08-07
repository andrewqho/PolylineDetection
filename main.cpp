#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/flann.hpp>
#include <chrono>

/*
 main.cpp
 polygon_detection

 Created by Andrew Ho on 6/24/18.
 Copyright © 2018 Andrew Ho, California Institute of Technology. All rights reserved.
 */

#include "line_detection.hpp"

using namespace std::chrono;

int main(int argc, char* argv[]) {
    String sequence_name = "without_bg";
    String sequencepath = "/Users/aqho/Desktop/Pose_Aquisition/Polygon_Detection/Polygon_Detection/static/" + sequence_name + "/*.png";
    String outputpath = "/Users/aqho/Desktop/Pose_Aquisition/Polygon_Detection/Polygon_Detection/video_output/"+ sequence_name + "/";
    vector<String> filenames;
    cv::glob(sequencepath, filenames);

    // Drawing parameters
    int radius = 4;
    
    // Running parameters
    bool specific = false;
    size_t limiter = 200;
    
    if(specific) {
        limiter = 1;
    }

    if(filenames.size() < limiter){
        limiter = filenames.size();
    }

    // Tuning parameters
    bool backgroundPresent = false;
    bool croppingNeeded = false;

    cout << sequence_name << endl;

    for (size_t image_num = 0; image_num<limiter; image_num++){
        
        Mat src;

        if(!specific) {
            src = imread(filenames[image_num]);
        }
        
        else {
            src = imread("/Users/aqho/Desktop/Pose_Aquisition/Polygon_Detection/Polygon_Detection/static/without_bg/damico.png", 1);
        }

        // Crop for specific images that have different text/dimensions
        if(croppingNeeded){
            cv::Rect init_ROI(0, 100, 1200, 550);
            src = src(init_ROI);
        }

        // =============================================================================
        // FIND AND FIX BOUNDING BOX

        // Blur image
        Mat src_blur;
        GaussianBlur(src, src_blur, Size(5, 5), 0, 0);
        
        Rect ROI = findROI(src_blur, backgroundPresent);
        // =============================================================================
        // Crop image

        Mat cropped_WGE = cropImg(src, ROI, backgroundPresent);

        // =============================================================================
        // Prepare inputs

        // Prepare input for ST
        Mat cropped_src = crop(src, ROI);

        // Resize images
        Mat resized_WGE;
        cv::resize(cropped_WGE, resized_WGE, cv::Size(800, 800));
        
        Mat resized_src;
        cv::resize(cropped_src, resized_src, cv::Size(800, 800));
        
        Mat resized_prewitt = prewittOperator(resized_src, 1.0);
        
        cvtColor(resized_src, resized_src, COLOR_BGR2GRAY);
        
        // =============================================================================
        // THRESHOLD

        // Calculate initial Otsu's Value
        int otsuVal = getOtsuVal(resized_prewitt);
        
        // High pass threshold
        int inc_otsuVal = otsuVal + 0.2*(255 - otsuVal);
        Mat thresh = global_thresh(resized_prewitt, inc_otsuVal, 0);
        
        cvtColor(thresh, thresh, COLOR_BGR2GRAY);
        
        morphOp(&thresh, 0, 1);
        morphOp(&thresh, 1, 0);
        
        // =============================================================================
        // HOUGH LINE TRANSFORM
        vector<Vec4f> hough_lines;
        
        HoughLinesP(thresh, hough_lines, 1, CV_PI/180, 50, 20, 10 );
        
        Mat hough_dst = Mat::zeros(resized_src.size(), CV_8UC3 );
        for( size_t i = 0; i < hough_lines.size(); i++ )
        {
            Vec4i l = hough_lines[i];
            line( hough_dst, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0,0,255), 3);
        }

        // =============================================================================
        // LINE SEGMENT DETECTION
        cv::Ptr<cv::LineSegmentDetector> det;
        det = cv::createLineSegmentDetector();

        vector<Vec4f> LSD_lines;

        det->detect(thresh, LSD_lines);

        Mat LSD_dst = Mat::zeros(resized_src.size(), CV_8UC3 );
        det->drawSegments(LSD_dst, LSD_lines);
        // =============================================================================
        // CORNER DETECTION
        
        vector<Point2f> ST_points = goodFeaturesToTrack_Callback(30, thresh);

        Mat ST_dst = resized_src.clone();
        for( size_t i = 0; i < ST_points.size(); i++ )
        {
            circle(ST_dst, ST_points[i], radius, Scalar(0, 0, 255), FILLED );
        }
        
        
        // =============================================================================
        // FILTER POINTS

        cv::Mat_<float> Hough_points = convertPoints(hough_lines);
        cv::Mat_<float> LSD_points = convertPoints(LSD_lines);

        cv::flann::Index Hough_tree = createKDTree(Hough_points);
        cv::flann::Index LSD_tree = createKDTree(LSD_points);

        // Hold filtered_corners
        vector<gNode> temp_filtered_corners = filterPoints(ST_points, Hough_points, &Hough_tree, LSD_points, &LSD_tree);

        // Convert to pointers
        vector<gNode*> filtered_corners;
        for( size_t i = 0; i < temp_filtered_corners.size(); i++ )
        {
            filtered_corners.push_back(&temp_filtered_corners[i]);
        }

        // Draw filtered points
        
        Mat final_img = resized_src.clone();
        cvtColor(final_img, final_img, COLOR_GRAY2BGR);
        for(auto && tpoint: filtered_corners)
        {
            circle(final_img, Point(tpoint->x, tpoint->y), radius, Scalar(0, 0, 255), -1);
        }

        string final_filename = outputpath + to_string(image_num) + "_final.png";
        imwrite(final_filename, final_img);

        imshow("Filtered Points", final_img);

        // =============================================================================
        // FIND ALL POSSIBLE EDGES

        // Prepare Frame
        
        // High pass threshold
        Mat frame = global_thresh(resized_prewitt, inc_otsuVal, 0);
        
        cvtColor(frame, frame, COLOR_BGR2GRAY);
        
        morphOp(&frame, 0, 1);
        morphOp(&frame, 1, 0);
        
        imshow("Frame", frame);
        
        // Find edges
        findEdges(frame, filtered_corners);
        
        Mat filtered = resized_src.clone();
        cvtColor(filtered, filtered, COLOR_GRAY2BGR);
        for(auto && tpoint: filtered_corners)
        {
            for(auto && neighbor: tpoint->adjacent)
            {
                line(filtered, Point(tpoint->x, tpoint->y), Point(neighbor->x, neighbor->y), Scalar(0, 255, 0), 2);
            }
        }
        
        imshow("Constructed Edges", filtered);
        
        // =============================================================================
        // FILTER DETECTED EDGES
        
//        filterEdges(filtered_corners);

        restrictEdges(filtered_corners);
        filterEdges(filtered_corners);
        filterEdges(filtered_corners);
        
        Mat filtered_final = resized_src.clone();
        
        cvtColor(filtered_final, filtered_final, COLOR_GRAY2BGR);
        for(auto && tpoint: filtered_corners)
        {
            circle(filtered_final, Point(tpoint->x, tpoint->y), radius, Scalar(255, 0, 255), -1);
        }
        for(auto && tpoint: filtered_corners)
        {
            for(auto && neighbor: tpoint->adjacent)
            {
                line(filtered_final, Point(tpoint->x, tpoint->y), Point(neighbor->x, neighbor->y), Scalar(255, 255, 0), 2);
            }
        }
        imshow("Filtered edges", filtered_final);
        waitKey(0);
    }
}
