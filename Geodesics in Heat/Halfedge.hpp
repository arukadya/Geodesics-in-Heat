//
//  Halfedge.hpp
//  m1morimori
//
//  Created by sunouchitoshiki on 2023/04/25.
//

#ifndef Halfedge_hpp
#define Halfedge_hpp
#include "Eigen/Core"
#include "Eigen/LU"
#include <filesystem>
#include <set>
#include <map>
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
struct Halfedge{
    int index;
    int face;
    int h_next;
    int h_prev;
    int v_src;
    int v_tgt;
    Halfedge(int i,int src,int tgt);
    Halfedge();
    int h_opp();
    void setNextPrev(int next,int prev);
    void setFace(int f);
};
struct Mesh{
    std::vector<Eigen::Vector3d> Vertices;
    std::vector<std::vector<int>> Faces;
    std::vector<Halfedge>HalfedgeList;
    std::vector<int>v2he;
    Mesh(std::vector<Eigen::Vector3d> &V,
         std::vector<std::vector<int>> &F);
    Mesh();
    void outputOFF(const char* OutputFileName);
    void outputOBJ(const char* OutputFileName);
    void input(std::string InputFlieName);
    int inputOBJ(std::string InputFlieName);
    void inputOFF(const char* InputFileName);
    void makeHalfedgeList();
    void printOneLink(int v_index);
    int h_cw(int h_index);
    int h_ccw(int h_index);
    double angle_deficit(int v_index);
    double sum_angle_deficit();
    double ave_edge_length();
};

#endif /* Halfedge_hpp */
