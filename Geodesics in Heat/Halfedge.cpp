//
//  Halfedge.cpp
//  m1morimori
//
//  Created by sunouchitoshiki on 2023/04/25.
//

#include "Halfedge.hpp"
Halfedge::Halfedge(int i,int src,int tgt){
    index = i;
    face = -1;
    v_src = src;
    v_tgt = tgt;
    h_next = -1;
    h_prev = -1;
}
Halfedge::Halfedge(){
    index = -1;
    face = -1;
    v_src = -1;
    v_tgt = -1;
    h_next = -1;
    h_prev = -1;
}
int Halfedge::h_opp(){
    if(index %2 == 0)return index+1;
    else return index-1;
}
void Halfedge::setNextPrev(int next,int prev){
    h_next = next;
    h_prev = prev;
}
void Halfedge::setFace(int f){
    face = f;
}
int Mesh::h_cw(int h_index){
    return HalfedgeList[HalfedgeList[h_index].h_opp()].h_next;
}
int Mesh::h_ccw(int h_index){
    if(HalfedgeList[h_index].h_prev == -1){
        //std::cout << "error_ccw" << std::endl;
        return -1;
    }
    return HalfedgeList[HalfedgeList[h_index].h_prev].h_opp();
}

void Mesh::input(std::string InputFileName){
    std::cout << "InputFileName:" << InputFileName << std::endl;
    if(std::filesystem::path(InputFileName).extension() == ".obj")inputOBJ(InputFileName);
    if(std::filesystem::path(InputFileName).extension() == ".off")inputOFF(InputFileName.c_str());
}
int Mesh::inputOBJ(std::string InputFileName){
    std::ifstream Inputfile(InputFileName);
    if (!Inputfile.is_open()) {
        std::cerr << "Could not open the file - '"
             << InputFileName << "'" << std::endl;
        return EXIT_FAILURE;
    }
    std::string line;
    std::string word;
    
    while(std::getline(Inputfile,line)){
        std::stringstream ss_nums{line};
        getline(ss_nums,word,' ');
        if(word[0] == '#' || word[0] == '\n')continue;
        if(word[0] == 'v' && word.size() == 1){
            //std::cout << word << ",";
            Eigen::Vector3d v0;
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            //std::cout << word << ",";
            v0.x() = std::stof(word);
            
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            v0.y() = std::stof(word);
            //std::cout << word << ",";
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            v0.z() = std::stof(word);
            //std::cout << word << ",";
            Vertices.push_back(v0);
            //std::cout << std::endl;
        }
        if(word[0] == 'f'){
            std::vector<int> f;
            while(getline(ss_nums,word,' ')){
                f.push_back(std::stoi(word)-1);
            };
            Faces.push_back(f);
        }
    };
    Inputfile.close();
    return EXIT_SUCCESS;
}
void Mesh::inputOFF(const char* InputFileName){
    FILE *ifp = fopen(InputFileName,"r");
    int num_vertices, num_faces,dummy;
    fscanf(ifp, "OFF %d %d %d", &num_vertices, &num_faces, &dummy);
    for(int i=0;i<num_vertices;i++){
        double x,y,z;//点の入力
        fscanf(ifp, "%lf %lf %lf", &x, &y, &z);
        Eigen::Vector3d v(x,y,z);
        Vertices.push_back(v);
    }
    for(int i=0;i<num_faces;i++){//面の入力
        int num_size, v0, v1, v2;
        fscanf(ifp, "%d %d %d %d", &num_size, &v0, &v1, &v2);
        Faces.push_back({v0,v1,v2});
    }
    fclose(ifp);
}
void Mesh::outputOFF(const char* OutputFileName){
    FILE *ofp = fopen(OutputFileName,"w");
    fprintf(ofp,"OFF\n");
    fprintf(ofp,"%lu %lu 0\n",Vertices.size(),Faces.size());
    for(auto &x:Vertices){
        fprintf(ofp,"%lf %lf %lf\n",x(0),x(1),x(2));
    }
    for(int i=0;i<Faces.size();i++){
        fprintf(ofp,"3 ");
        for(int j=0;j < 3;j++){
            fprintf(ofp,"%d",Faces[i][j]);
            if(j != 2)fprintf(ofp," ");
            else fprintf(ofp,"\n");
        }
    }
    fclose(ofp);
}
void Mesh::outputOBJ(const char* OutputFileName){
    FILE *ofp = fopen(OutputFileName,"w");
    for(auto &x:Vertices){
        fprintf(ofp,"v %lf %lf %lf\n",x(0),x(1),x(2));
    }
    for(int i=0;i<Faces.size();i++){
        fprintf(ofp,"f ");
        for(int j=0;j < 3;j++){
            fprintf(ofp,"%d",Faces[i][j]+1);
            if(j != 2)fprintf(ofp," ");
            else fprintf(ofp,"\n");
        }
    }
    fclose(ofp);
}
Mesh::Mesh(std::vector<Eigen::Vector3d> &V,
           std::vector<std::vector<int>> &F){
    Vertices = V;
    Faces = F;
}
Mesh::Mesh(){
}
void Mesh::makeHalfedgeList(){
    std::map<std::pair<int, int>,int>map;
    int halfedgecnt = 0;
    for(int i=0;i<Faces.size();i++){
        for(int j=0;j<Faces[i].size();j++){
            int v0 = Faces[i][j%Faces[i].size()];
            int v1 = Faces[i][(j+1)%Faces[i].size()];
            std::pair<int, int>key_now = std::make_pair(v0, v1);
            std::pair<int, int>key_now_opp = std::make_pair(v1, v0);
            int h_now_index = halfedgecnt;
            if(map.count(key_now) == 0 && map.count(key_now_opp) == 0){
                map.emplace(key_now,halfedgecnt);
                map.emplace(key_now_opp,halfedgecnt+1);
                Halfedge h_now = Halfedge(h_now_index,v0,v1);
                Halfedge h_opp = Halfedge(h_now_index+1,v1,v0);
                h_now.setFace(i);
                HalfedgeList.push_back(h_now);
                HalfedgeList.push_back(h_opp);
                halfedgecnt+=2;
            }
            else HalfedgeList[map.at(key_now)].setFace(i);
        }
    }
    for(int i=0;i<Faces.size();i++){
        for(int j=0;j<Faces[i].size();j++){
            int v0 = Faces[i][j%Faces[i].size()];
            int v1 = Faces[i][(j+1)%Faces[i].size()];
            int v2 = Faces[i][(j+2)%Faces[i].size()];
            int v3 = Faces[i][(j+3)%Faces[i].size()];
            std::pair<int, int>key_prev = std::make_pair(v0, v1);
            std::pair<int, int>key_now = std::make_pair(v1, v2);
            std::pair<int, int>key_next = std::make_pair(v2, v3);
            int h_now_index = map[key_now];
            int h_next_index = map[key_next];
            int h_prev_index = map.at(key_prev);
            HalfedgeList[h_now_index].setNextPrev(h_next_index, h_prev_index);
        }
    }
    v2he.resize(Vertices.size());
    for(auto &h:HalfedgeList){
        v2he[h.v_src] = h.index;
    }
}
void Mesh::printOneLink(int v_index){
    Halfedge h_now = HalfedgeList[v2he[v_index]];
    int start = h_now.v_tgt;
    do{
        if(h_cw(h_now.index) == -1)break;
        h_now = HalfedgeList[h_cw(h_now.index)];
    }while(h_cw(h_now.index) == -1 && h_now.v_tgt != start);
    std::cout << std::endl;
    start = h_now.v_tgt;
    do{
        if(h_ccw(h_now.index) == -1)break;
        h_now = HalfedgeList[h_ccw(h_now.index)];
    }while(h_ccw(h_now.index) != -1 && h_now.v_tgt != start);
}
double Mesh::angle_deficit(int v_index){
    double sum = 0;
    Halfedge h_now = HalfedgeList[v2he[v_index]];
    int start = h_now.v_tgt;
//    std::cout << std::endl << "start_h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" << h_now.v_tgt <<std::endl;
    Eigen::Vector3d x0 = Vertices[h_now.v_src];
    do{
        if(h_cw(h_now.index) == -1)break;
        h_now = HalfedgeList[h_cw(h_now.index)];
//        std::cout << "h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" <<   h_now.v_tgt <<std::endl;
    }while(h_cw(h_now.index) != -1 && h_now.v_tgt != start);
//    std::cout << "fin_cw" << std::endl << std::endl;
    start = h_now.v_tgt;
    do{
        if(h_ccw(h_now.index) == -1)break;
        Eigen::Vector3d x1 = Vertices[h_now.v_tgt] - x0;
//        std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
        h_now = HalfedgeList[h_ccw(h_now.index)];
        Eigen::Vector3d x2 = Vertices[h_now.v_tgt] - x0;
        double cos = x1.dot(x2)/(x1.norm()*x2.norm());
        sum+= acos(cos);
    }while(h_ccw(h_now.index) != -1 && h_now.v_tgt != start);
//    std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
//    std::cout << "fin_ccw" << std::endl << std::endl;
    return 2*M_PI - sum;
}
double Mesh::sum_angle_deficit(){
    double sum = 0;
    for(int i=0;i<Vertices.size();++i){
        sum += angle_deficit(i);
    }
    return sum;
}
double Mesh::ave_edge_length(){
    double sum_not_bound = 0;
    double sum_bound = 0;
    for(int i=0;i<HalfedgeList.size();++i){
        Eigen::Vector3d v0 = Vertices[HalfedgeList[i].v_src];
        Eigen::Vector3d v1 = Vertices[HalfedgeList[i].v_tgt];
        if(HalfedgeList[i].h_next != -1)sum_not_bound += (v0 - v1).norm();
        else sum_bound += (v0 - v1).norm();
    }
    return sum_not_bound/2 + sum_bound;
}
