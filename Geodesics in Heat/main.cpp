//
//  main.cpp
//  Geodesics in Heat
//
//  Created by 須之内俊樹 on 2023/07/20.
//

#include <iostream>
#include "Halfedge.hpp"
int main(int argc, const char * argv[]) {
    Mesh test;
    if(argc == 1)test.input("bun_zipper.off");
    else if(argc >= 2)test.input(argv[1]);
    std::cout << test.Vertices.size() << "," << test.Faces.size() << std::endl;
    return 0;
}
