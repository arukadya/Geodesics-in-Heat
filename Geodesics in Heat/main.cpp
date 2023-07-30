//
//  main.cpp
//  Geodesics in Heat
//
//  Created by 須之内俊樹 on 2023/07/20.
//

#include <iostream>
#include "Laplacian_Mesh.hpp"
int main(int argc, const char * argv[]) {
    Laplacian_Mesh test;
    std::vector<int>delta_id = {0};
    std::string inputFileName;
    std::string outputFileName;
    if(argc == 1){
        inputFileName = "294_kitten_uniform.off";
        outputFileName = "kitten";
    }
    else if(argc >= 2){
        inputFileName = argv[1];
        outputFileName = "test";
    }
    test.input(inputFileName);
    test.makeHalfedgeList();
    double t = test.ave_edge_length();
    std::cout << "t=" << t << std::endl;
    test.cal_geodescis_distance(t*t, delta_id);
    test.outputVTK((outputFileName + "_geodescis_distance.vtk").c_str());
    test.output_Vector_VTK((outputFileName + "_gradient.vtk").c_str());
    return 0;
}
