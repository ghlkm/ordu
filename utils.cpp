//
// Created by 12859 on 2023/5/22.
//

#include "utils.h"
#include <iostream>
using namespace std;

void kskyband_read(const string &filename, vector<vector<int>> &ret){
    FILE *in=fopen(filename.c_str(), "r");
    int cur;
    char UNUSED;
    int flag=1;
    vector<int> one_onion;
    while(flag){
        flag=fscanf(in, "%d", &cur);
        if(feof(in)){
            break;
        }
        if(flag){
            one_onion.push_back(cur);
        }else{
            flag=fscanf(in, "%c", &UNUSED);
            if(flag){
                if(UNUSED=='#'){
                    flag=fscanf(in, "%d", &cur);
                }
                flag=fscanf(in, "%c", &UNUSED);
                if(cur!=1){
                    ret.push_back(one_onion);
                    one_onion=vector<int>();
                }
            }// else eof
        }
    }
    fclose(in);
    ret.push_back(one_onion);
    cout<<ret.size()<<endl;
    for(auto &i:ret){
        cout<<i.size()<<endl;
    }
}

