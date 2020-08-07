/*
 line_detection.hpp
 polygon_detection
 
 Created by Andrew Ho on 6/24/18.
 Copyright © 2018 Andrew Ho. All rights reserved.
 */

#ifndef line_detection_hpp
#define line_detection_hpp

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <functional>

#include <opencv2/line_descriptor.hpp>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"
#include <opencv2/flann.hpp>

using namespace cv;
using namespace std;

struct gNode{
    // Positional Info
    float x;
    float y;
    
    // Polygon Data
    
    // Hold pointers to all adjacent nodes
    vector<gNode*> adjacent;
    
    // Denote the chain the gNode belongs to
    int chain_index;
    
    // Mark if this point has been rejected before
    bool rejected;
    
    gNode(Point2f in_pt){
        // Initialize positional information
        x = in_pt.x;
        y = in_pt.y;
        
        // Initialize chain_index to -1
        chain_index = -1;
        
        // Initialize rejected status to false
        rejected = false;
        
    }
    
    gNode(float in_x, float in_y){
        // Initialize positional information
        x = in_x;
        y = in_y;
        
        // Initialize chain_index to -1
        chain_index = -1;
        
        // Initialize rejected status to false
        rejected = false;
    }
    
    ~gNode(){
        
    }
    
    void addEdge(gNode *vertex){
        adjacent.push_back(vertex);
    }
    
    void removeEdge(gNode *vertex){
        for(int i = 0; i < adjacent.size(); i++){
            if(adjacent[i] == vertex){
                adjacent.erase(adjacent.begin() + i);
                break;
            }
        }
    }
    
    bool edgeExist(gNode *vertex){
        for(int i = 0; i < adjacent.size(); i++){
            if(adjacent[i] == vertex){
                return true;
            }
        }
        return false;
    }
    
    int getNumNeighbors(){
        return adjacent.size();
    }
    
    void printNeighbors(){
        for(auto && neighbors: adjacent){
            cout << "(" << neighbors->x << ", " << neighbors->y << ")";
        }
    }
    
    void reset(){
        chain_index = -1;
        for(auto && adj_node : adjacent){
            adj_node->removeEdge(this);
        }
        adjacent.clear();
    }
    
    String printCoords(){
        String ret_str  = "(" + to_string(x) + ", " + to_string(y) + ") ";
        return ret_str;
    }
};

double getOtsuVal( const Mat& );
Mat global_thresh(Mat, int, int);
Mat thresh_callback(int, void*, Mat);
float calcMagnitude(Point, Point) ;
float dotproduct(vector<float>, vector<float>);
float angleBetween(Point, Point, Point);
vector<Point2f> goodFeaturesToTrack_Callback(int, Mat);
tuple<float, vector<Point2f>> calcCost(vector<Point2f>, Point2f);
vector<Point2f> TSA(vector<Point2f>);
Mat adaptive_thresh(Mat, int, int, int);
void thinning(Mat& im);
void thinningIteration(Mat& im, int);
Point2f calcCenter(Point2f, Point2f, Point2f);
Mat prewittOperator(Mat, float);
Mat WGE(Mat, float);
Mat cropImg(Mat, Rect, bool);
Mat crop(Mat, Rect);
Rect findROI(Mat, bool);
cv::Mat_<float> convertPoints(vector<Vec4f> input_lines);
cv::flann::Index createKDTree(cv::Mat_<float> flannPoints);
vector<gNode> filterPoints(vector<Point2f>, cv::Mat_<float>, cv::flann::Index*, cv::Mat_<float>, cv::flann::Index*);
Vec2f calcPixGradient(Mat, int, int);
void findEdges(Mat, vector<gNode*>);
void morphOp(Mat*, int, int);
void filterEdges(vector<gNode*>);
float calcMagnitudeGNode(gNode*, gNode*);
void restrictEdges(vector<gNode*>);
bool testPixels(Mat, int, int, int);
#endif /* line_detection_hpp */

