//
// Created by 12859 on 2023/5/22.
//

#ifndef IPREF_UTILS_H
#define IPREF_UTILS_H

#include<string>
#include<vector>
void kskyband_read(const std::string &filename, std::vector<std::vector<int>> &ret);

// inside [0, 1]
template<typename T>
inline bool legal(const std::vector<T> &v){
    for (const T& val: v) {
        if(val >1 || val <0 ) return false;
    }
    return true;
}




#endif //IPREF_UTILS_H
