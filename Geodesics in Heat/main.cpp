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
    if(argc == 1)test.input("bun_zipper.off");
    else if(argc >= 2)test.input(argv[1]);
    test.makeHalfedgeList();
    test.cal_TriArea();
    test.cal_Laplacian();
//    std::cout << test.Laplacian << std::endl;
//    std::cout << test.Laplacian_C << std::endl;
//    for(int i=0;i<test.vertex_TriArea.size();++i)std::cout << i << "," << test.vertex_TriArea[i] << std::endl;
    return 0;
}
