//  Created by sunouchitoshiki on 2023/04/25.
//

#include "Laplacian_Mesh.hpp"
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
int Laplacian_Mesh::h_cw(int h_index){
    return HalfedgeList[HalfedgeList[h_index].h_opp()].h_next;
}
int Laplacian_Mesh::h_ccw(int h_index){
    if(HalfedgeList[h_index].h_prev == -1){
        //std::cout << "error_ccw" << std::endl;
        return -1;
    }
    return HalfedgeList[HalfedgeList[h_index].h_prev].h_opp();
}

void Laplacian_Mesh::input(std::string InputFileName){
    std::cout << "InputFileName:" << InputFileName << std::endl;
    if(std::filesystem::path(InputFileName).extension() == ".obj")inputOBJ(InputFileName);
    if(std::filesystem::path(InputFileName).extension() == ".off")inputOFF(InputFileName.c_str());
}
int Laplacian_Mesh::inputOBJ(std::string InputFileName){
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
void Laplacian_Mesh::inputOFF(const char* InputFileName){
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
void Laplacian_Mesh::outputOFF(const char* OutputFileName){
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
void Laplacian_Mesh::outputOBJ(const char* OutputFileName){
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
void Laplacian_Mesh::outputVTK(const char* OutputFileName){
    FILE *ofp = fopen(OutputFileName,"w");
    fprintf(ofp, "# vtk DataFile Version 2.0\n");
    fprintf(ofp, "Title Data\n");
    fprintf(ofp, "ASCII\n");
    fprintf(ofp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(ofp, "POINTS           %d float\n",Vertices.size());
    for(auto &x:Vertices)fprintf(ofp,"%lf %lf %lf\n",x(0),x(1),x(2));
    fprintf(ofp, "CELLS %d %d\n", Faces.size() ,Faces.size()*4);
    for(int i=0;i<Faces.size();i++){
        fprintf(ofp,"3 ");
        for(int j=0;j < 3;j++){
            fprintf(ofp,"%d",Faces[i][j]);
            if(j != 2)fprintf(ofp," ");
            else fprintf(ofp,"\n");
        }
    }
    fprintf(ofp,"CELL_TYPES   %d\n",Faces.size());
    for(int i=0;i<Faces.size();i++)fprintf(ofp,"5\n");
    fprintf(ofp,"POINT_DATA   %d\n",geodescis_distance.size());
    fprintf(ofp, "SCALARS Temp float\n");
    fprintf(ofp, "LOOKUP_TABLE default\n");
    for(auto &x:geodescis_distance)fprintf(ofp, "%lf\n",x);
    fclose(ofp);
}
void Laplacian_Mesh::output_Vector_VTK(const char* OutputFileName){
    FILE *ofp = fopen(OutputFileName,"w");
    fprintf(ofp, "Title Data\n");
    fprintf(ofp, "ASCII\n");
    fprintf(ofp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(ofp, "POINTS           %d float\n",Faces.size());
    
    for(int i=0;i<Faces.size();++i){
        Eigen::Vector3d center_Face = (Vertices[Faces[i][0]] + Vertices[Faces[i][1]] + Vertices[Faces[i][2]])/3;
        fprintf(ofp,"%lf %lf %lf\n",center_Face.x(),center_Face.y(),center_Face.z());
    }
    fprintf(ofp, "CELLS %d %d\n", Faces.size() ,Faces.size());
    for(int i=0;i<Faces.size();++i){
        fprintf(ofp,"1 %d",i);
        fprintf(ofp,"\n");
    }
//    fprintf(ofp,"CELL_TYPES   %d\n",Faces.size());
    
    for(int i=0;i<Faces.size();i++)fprintf(ofp,"1\n");
    fprintf(ofp, "POINT_DATA    %d\n",Faces.size());
    fprintf(ofp, "Vector grad float\n");
    for(int i=0;i<gradients.size();++i){
        Eigen::Vector3d u = gradients[i];
        fprintf(ofp, "%d %lf %lf %lf\n",i,u.x(),u.y(),u.z());
    }
    fclose(ofp);
}
Laplacian_Mesh::Laplacian_Mesh(std::vector<Eigen::Vector3d> &V,
           std::vector<std::vector<int>> &F){
    Vertices = V;
    Faces = F;
}
Laplacian_Mesh::Laplacian_Mesh(){
}
void Laplacian_Mesh::makeHalfedgeList(){
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
void Laplacian_Mesh::printOneLink(int v_index){
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
double Laplacian_Mesh::angle_deficit(int v_index){
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
double Laplacian_Mesh::sum_angle_deficit(){
    double sum = 0;
    for(int i=0;i<Vertices.size();++i){
        sum += angle_deficit(i);
    }
    return sum;
}
void Laplacian_Mesh::cal_TriArea(){
//    std::cout << "cal_TriArea" << std::endl;
    vertex_TriArea.resize(Vertices.size());
    for(int i=0;i<Vertices.size();++i){
        double sum = 0;
        Halfedge h_now = HalfedgeList[v2he[i]];
        int start = h_now.v_tgt;
//        std::cout << std::endl << "start_h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" << h_now.v_tgt <<std::endl;
        Eigen::Vector3d x0 = Vertices[h_now.v_src];
        do{
            if(h_cw(h_now.index) == -1)break;
            h_now = HalfedgeList[h_cw(h_now.index)];
//            std::cout << "h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" <<   h_now.v_tgt <<std::endl;
        }while(h_cw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "fin_cw" << std::endl << std::endl;
        start = h_now.v_tgt;
        int cnt = 0;
        do{
            if(h_ccw(h_now.index) == -1)break;
            Eigen::Vector3d x1 = Vertices[h_now.v_tgt] - x0;
//            std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
            h_now = HalfedgeList[h_ccw(h_now.index)];
            Eigen::Vector3d x2 = Vertices[h_now.v_tgt] - x0;
            double area = x1.cross(x2).norm()/6;
            sum += area;
        }while(h_ccw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
//        std::cout << "fin_ccw" << std::endl << std::endl;
        vertex_TriArea[i] = sum;
//        std::cout << "area_" << i << "=" << sum << std::endl;
    }
}
void Laplacian_Mesh::cal_Laplacian(){
    cal_TriArea();
    Laplacian.resize(Vertices.size(), Vertices.size());
    Laplacian_C.resize(Vertices.size(), Vertices.size());
    Laplacian_C_heat.resize(Vertices.size(), Vertices.size());
    std::vector<Triplet> triplets;
    std::vector<Triplet> triplets_C;
    std::vector<Triplet> triplets_A;
    for(int i=0;i<Vertices.size();++i){
        double sum = 0;
        double sum_area = 0;
        Halfedge h_now = HalfedgeList[v2he[i]];
        int start = h_now.v_tgt;
//        std::cout << std::endl << "start_h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" << h_now.v_tgt <<std::endl;
        Eigen::Vector3d x0 = Vertices[h_now.v_src];
        do{
            if(h_cw(h_now.index) == -1)break;
            h_now = HalfedgeList[h_cw(h_now.index)];
//            std::cout << "h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" <<   h_now.v_tgt <<std::endl;
        }while(h_cw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "fin_cw" << std::endl << std::endl;
        start = h_now.v_tgt;
        do{
            int j = h_now.v_tgt;
            if(h_ccw(h_now.index) == -1)break;
            Eigen::Vector3d a = Vertices[h_now.v_tgt];
            Eigen::Vector3d b = Vertices[HalfedgeList[h_now.h_next].v_tgt];
            Halfedge h_opp = HalfedgeList[h_now.h_opp()];
            Eigen::Vector3d c = Vertices[HalfedgeList[h_opp.h_next].v_tgt];
            Eigen::Vector3d x1 = x0 - b;
            Eigen::Vector3d x2 = a - b;
            Eigen::Vector3d x3 = x0 - c;
            Eigen::Vector3d x4 = a - c;
            double tan_alpha = x1.cross(x2).norm() / x1.dot(x2);
            double tan_beta = x3.cross(x4).norm() / x3.dot(x4);
//            double w_ij = (1/tan_alpha + 1/tan_beta)/vertex_TriArea[i];
            double w_ij = (1.0/tan_alpha + 1.0/tan_beta);
            sum += w_ij;
            sum_area += vertex_TriArea[i];
            triplets.emplace_back(i,j,w_ij);
//            triplets_C.emplace_back(i,j,w_ij*vertex_TriArea[i]);
            triplets_C.emplace_back(i,j,w_ij/2);
            h_now = HalfedgeList[h_ccw(h_now.index)];
        }while(h_ccw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
//        std::cout << "fin_ccw" << std::endl << std::endl;
        triplets.emplace_back(i,i, -sum);
//        triplets_C.emplace_back(i,i, -sum*sum_area);
        triplets_C.emplace_back(i,i, -sum/2);
        triplets_A.emplace_back(i,i, sum_area);
    }
    Laplacian.setFromTriplets(triplets.begin(), triplets.end());
    Laplacian_C.setFromTriplets(triplets_C.begin(), triplets_C.end());
    Laplacian_C_heat.setFromTriplets(triplets_A.begin(), triplets_A.end());
}
void Laplacian_Mesh::set_deltaOne(std::vector<int> &vertex_ids){
//    std::cout << Vertices.size() << std::endl;
    delta.resize(Vertices.size());
    for(auto &x:delta)x = 0;
//    for(auto &x:vertex_ids){
//        std::cout << x << std::endl;
//        delta[x] = 1;
//        std::cout << x << std::endl;
//    }
    for(int i=0;i<vertex_ids.size();++i){
        delta(vertex_ids[i]) = 1;
    }
}
void Laplacian_Mesh::cal_heat(double t){
//    std::cout << Laplacian_C << std::endl;
//    heat.resize(Vertices.size());
//    std::cout << Laplacian_C_heat << std::endl;
    Laplacian_C_heat -= t*Laplacian_C;
//    std::cout << Laplacian_C_heat << std::endl;
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(Laplacian_C_heat);
    heat = solver.solve(delta);
//    std::cout << heat.transpose() << std::endl;
}

void Laplacian_Mesh::cal_gradient(){
    gradients.resize(Faces.size());
    for(int i=0;i<Faces.size();++i){
        gradients[i] = {0,0,0};
        for(int j=0;j<3;++j){
            Eigen::Vector3d a = Vertices[Faces[i][j]];
            Eigen::Vector3d b = Vertices[Faces[i][(j+1)%3]];
            Eigen::Vector3d c = Vertices[Faces[i][(j+2)%3]];
            Eigen::Vector3d e = c-b;
            Eigen::Vector3d Normal = (b-a).cross(c-a);
            double area_f = Normal.norm()/2;
            gradients[i] += 1.0/(2*area_f)*heat[Faces[i][j]]*(Normal.normalized()).cross(e);
        }
    }
}
void Laplacian_Mesh::cal_integrated_div(){
    integrated_div.resize(Vertices.size());
    for(int i=0;i<Vertices.size();++i){
        integrated_div[i] = 0;
        Halfedge h_now = HalfedgeList[v2he[i]];
        int start = h_now.v_tgt;
//        std::cout << std::endl << "start_h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" << h_now.v_tgt <<std::endl;
        Eigen::Vector3d x0 = Vertices[h_now.v_src];
        do{
            if(h_cw(h_now.index) == -1)break;
            h_now = HalfedgeList[h_cw(h_now.index)];
//            std::cout << "h_index=" << h_now.index <<",h_src="<< h_now.v_src << ",h_tgt=" <<   h_now.v_tgt <<std::endl;
        }while(h_cw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "fin_cw" << std::endl << std::endl;
        start = h_now.v_tgt;
        do{
            int j = h_now.v_tgt;
            if(h_ccw(h_now.index) == -1)break;
            Eigen::Vector3d a = Vertices[h_now.v_tgt];
            Eigen::Vector3d b = Vertices[HalfedgeList[h_now.h_next].v_tgt];
            Eigen::Vector3d e1 = a-x0;
            Eigen::Vector3d e2 = b-x0;
            Eigen::Vector3d e3 = b-a;
            double tan1 = (-e2).cross(-e3).norm() / (-e2).dot(-e3);
            double tan2 = (-e1).cross(e3).norm() / (-e1).dot(e3);
//            std::cout << "div["<<i<<"] += " << "(1/2)*( ("<< 1.0/tan1 << ")*" << e1.dot(-gradients[h_now.face].normalized()) << " +  (" << 1.0/tan2 << ")*" << e2.dot(-gradients[h_now.face].normalized()) <<" ) = " << 0.5*( (1.0/tan1)*e1.dot(-gradients[h_now.face].normalized()) + (1.0/tan2)*e2.dot(-gradients[h_now.face].normalized()) ) << std::endl;
            integrated_div[i] += 0.5*( (1.0/tan1)*e1.dot(-gradients[h_now.face].normalized()) + (1.0/tan2)*e2.dot(-gradients[h_now.face].normalized()) );
            
            h_now = HalfedgeList[h_ccw(h_now.index)];
        }while(h_ccw(h_now.index) != -1 && h_now.v_tgt != start);
//        std::cout << "tgt=" << h_now.v_tgt << "," << "src=" << h_now.v_src << std::endl;
//        std::cout << "fin_ccw" << std::endl << std::endl;
    }
//    for(auto &x:integrated_div)std::cout << x << std::endl;
}
double Laplacian_Mesh::ave_edge_length(){
    double sum_not_bound = 0;
    double sum_bound = 0;
    int num_edges = Faces.size()/2*3;
    for(int i=0;i<HalfedgeList.size();++i){
        Eigen::Vector3d v0 = Vertices[HalfedgeList[i].v_src];
        Eigen::Vector3d v1 = Vertices[HalfedgeList[i].v_tgt];
        sum_not_bound += (v0 - v1).norm();
    }
    return sum_not_bound/num_edges/2;
}
void Laplacian_Mesh::cal_geodescis_distance(double t,std::vector<int> &vertex_ids){
    std::cout << "start" << std::endl;
    set_deltaOne(vertex_ids);
    std::cout << "set_delta" << std::endl;
    cal_TriArea();
    std::cout << "cal_TriArea" << std::endl;
//    for(auto &x:vertex_TriArea)std::cout << x << std::endl;
    cal_Laplacian();
    std::cout << "cal_Lapracian" << std::endl;
    cal_heat(t);
    std::cout << "cal_heat" << std::endl;
    cal_gradient();
    std::cout << "cal_gradient" << std::endl;
    cal_integrated_div();
    std::cout << "cal_integrated_div" << std::endl;
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(Laplacian_C);
//    Eigen::VectorXd r = Eigen::VectorXd::Random(Vertices.size());
//    std::cout << integrated_div.transpose() << std::endl;
    geodescis_distance = solver.solve(integrated_div);
//    for(auto &x:geodescis_distance)std::cout << x << std::endl;
//    std::cout << std::endl;
//    geodescis_distance = solver.solve(integrated_div.normalized());
//    for(auto &x:geodescis_distance)std::cout << x << std::endl;
//    geodescis_distance = solver.solve(r);
    double min = geodescis_distance.minCoeff();
    for(auto &x:geodescis_distance)x -= min;
//    for(auto &x:geodescis_distance)std::cout << x << std::endl;
}

