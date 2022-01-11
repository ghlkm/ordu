//
// Created by Keming Li on 2022/1/11.
//

#ifndef IPREF_REGION_H
#define IPREF_REGION_H
#include <vector>
using std::vector;
class region{
public:
    vector<int> topk;
    double radius;
    vector<vector<double>> cone;
public:
    region(vector<int> &topK, vector<vector<double>> &Cone){
        topk=topK;
        cone=Cone;
    }
    void setRadius(float new_r){
        this->radius=new_r;
    }
    void setRadius(double new_r){
        this->radius=new_r;
    }
};

#endif //IPREF_REGION_H
