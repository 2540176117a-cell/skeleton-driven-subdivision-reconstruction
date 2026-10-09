#include "StraightSkeletonSegmenter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <algorithm>
int main(int argc, char** argv) {
    if(argc!=3) { std::cerr<<"Usage: segment_demo skeleton.txt segments.obj\n"; return 1; }
    std::ifstream in(argv[1]); size_t nv=0, ne=0;
    if(!(in>>nv>>ne) || nv<2 || nv>1000000 || ne>10000000) return 2;
    std::vector<double> nodes(nv*3); std::vector<unsigned int> edges(ne*2);
    for(auto& v:nodes) if(!(in>>v) || !std::isfinite(v)) return 2;
    std::set<std::pair<unsigned int,unsigned int>> seen;
    for(size_t i=0;i<ne;++i) {
        auto& a=edges[2*i]; auto& b=edges[2*i+1];
        if(!(in>>a>>b) || a>=nv || b>=nv || a==b) return 2;
        auto key=std::minmax(a,b); if(!seen.emplace(key.first,key.second).second) return 2;
    }
    auto segments=StraightSkeletonSegmenter::ExtractStraightSegments(nodes,edges);
    std::ofstream out(argv[2]); if(!out) return 3;
    out<<std::setprecision(17); size_t base=1;
    for(const auto& s:segments) {
        out<<"# source_chain "<<s.sourceChainId<<"\n";
        for(size_t i=0;i<s.points.size();i+=3)
            out<<"v "<<s.points[i]<<" "<<s.points[i+1]<<" "<<s.points[i+2]<<"\n";
        out<<"l"; for(size_t i=0;i<s.points.size()/3;++i) out<<" "<<base+i;
        out<<"\n"; base+=s.points.size()/3;
    }
    if(!out) return 3;
    std::cout<<"segments="<<segments.size()<<"\n";
    return 0;
}
