//
//  Laplacian_Mesh.hpp
//  Geodesics in Heat
//
//  Created by 須之内俊樹 on 2023/07/20.
//

#ifndef Laplacian_Mesh_hpp
#define Laplacian_Mesh_hpp

#include "Eigen/Core"
#include "Eigen/SparseLU"
#include "Eigen/Dense"
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
#include <Eigen/IterativeLinearSolvers>

using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType,IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;

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
struct Laplacian_Mesh{
    std::vector<Eigen::Vector3d> Vertices;
    std::vector<Eigen::Vector3d> gradients;
    std::vector<std::vector<int>> Faces;
    std::vector<Halfedge>HalfedgeList;
    std::vector<int>v2he;
    Eigen::VectorXd heat;
    std::vector<double>geodescis_distance;
    Eigen::VectorXd delta;
    std::vector<double> vertex_TriArea;
    SparseMatrix Laplacian;
    SparseMatrix Laplacian_C;
    SparseMatrix Laplacian_C_heat;
    
    Laplacian_Mesh(std::vector<Eigen::Vector3d> &V,
         std::vector<std::vector<int>> &F);
    Laplacian_Mesh();
    void outputOFF(const char* OutputFileName);
    void outputOBJ(const char* OutputFileName);
    void outputVTK(const char* OutputFileName);
    void input(std::string InputFlieName);
    int inputOBJ(std::string InputFlieName);
    void inputOFF(const char* InputFileName);
    void makeHalfedgeList();
    void printOneLink(int v_index);
    int h_cw(int h_index);
    int h_ccw(int h_index);
    double angle_deficit(int v_index);
    double sum_angle_deficit();
    void cal_TriArea();
    void cal_Laplacian();
    void cal_heat(double t);
    void set_deltaOne(std::vector<int> vertex_id);
    void cal_gradient();
};




#endif /* Laplacian_Mesh_hpp */
