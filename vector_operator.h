
#ifndef UTK_BV_UTK_VECTOR_H
#define UTK_BV_UTK_VECTOR_H

#include <vector>
#include <cassert>
#include <cstddef>
#include <iostream>
using namespace std;

template<typename FLOAT1, typename FLOAT2>
vector<FLOAT1> operator+(const vector<FLOAT1> &v1, const vector<FLOAT2> &v2) {
    assert(v1.size() == v2.size());
    vector<FLOAT1> ret(v1);
    for (int i = 0; i < ret.size(); ++i) {
        ret[i] += v2[i];
    }
    return ret;
}

template<typename T>
vector<T> &operator+=(vector<T> &v1, const vector<T> &v2) {
    assert(v1.size() == v2.size());
    for (int i = 0; i < v1.size(); ++i) {
        v1[i] += v2[i];
    }
    return v1;
}

template<typename T>
vector<T> operator-(const vector<T> &v1, const vector<T> &v2) {
    assert(v1.size() == v2.size());
    vector<T> ret(v1);
    for (int i = 0; i < ret.size(); ++i) {
        ret[i] -= v2[i];
    }
    return ret;
}

template<typename T>
vector<T> &operator-=(vector<T> &v1, const vector<T> &v2) {
    assert(v1.size() == v2.size());
    for (int i = 0; i < v1.size(); ++i) {
        v1[i] -= v2[i];
    }
    return v1;
}

template<typename T>
T operator*(const vector<T> &v1, const vector<T> &v2) {
    assert(v1.size() == v2.size());
    T ret = T();
    for (int i = 0; i < v1.size(); ++i) {
        ret += v1[i] * v2[i];
    }
    return ret;
}

template<typename T>
vector<T> operator*(const T &num, const vector<T> &v) {
    vector<T> ret(v.size());
    for (int i = 0; i < v.size(); ++i) {
        ret[i] = num * v[i];
    }
    return ret;
}



template<typename T>
vector<T> &operator*=(vector<T> &v, const T &num) {
    for (int i = 0; i < v.size(); ++i) {
        v[i] *= v[i];
    }
    return v;
}

template<typename T>
vector<T> operator*(const vector<T> &v, const T &num) {
    return num * v;
}

template<typename T>
vector<T> operator/(const vector<T> &v, const T &num) {
    vector<T> ret(v.size());
    for (int i = 0; i < v.size(); ++i) {
        ret[i] = v[i] / num;
    }
    return ret;
}

template<typename T>
vector<T> &operator/=(vector<T> &v, const T &num) {
    for (int i = 0; i < v.size(); ++i) {
        v[i] = v[i] / num;
    }
    return v;
}

template<typename T>
ostream &operator<<(ostream &out, const vector<T> &v) {
    if (v.empty()) {
        return out;
    }
    out << v[0];
    for (auto i = 1; i < v.size(); ++i) {
        out << ", " << v[i];
    }
    return out;
}


#endif //UTK_BV_UTK_VECTOR_H
