
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

template<typename T1, typename T2>
vector<T1> operator-(const vector<T1> &v1, const vector<T2> &v2) {
    assert(v1.size() == v2.size());
    vector<T1> ret(v1);
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


template<typename T1, typename  T2>
T1 operator*(const vector<T1> &v1, const vector<T2> &v2) {
    assert(v1.size() == v2.size());
    T1 ret = T1();
    for (int i = 0; i < v1.size(); ++i) {
        ret += v1[i] * v2[i];
    }
    return ret;
}

template<typename T1, typename  T2>
T2 operator*(const T1 *v1, const vector<T2> &v2) {
    T2 ret = T2();
    for (int i = 0; i < v2.size(); ++i) {
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
vector<T> operator*(const int num, const vector<T> &v) {
    vector<T> ret(v.size());
    for (int i = 0; i < v.size(); ++i) {
        ret[i] = num * v[i];
    }
    return ret;
}

template<typename T>
vector<T> operator*(const T *num, const vector<T> &v) {
    T ret=0;
    for (int i = 0; i < v.size(); ++i) {
        ret += num[i] * v[i];
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

// Define a generic hash function for vector<T>
template<typename T>
struct VectorHash {
    size_t operator()(const vector<T>& vec) const {
        size_t seed = vec.size();
        for (const auto& elem : vec) {
            seed ^= hash<T>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template<typename T>
T abs(T val){
    if(val>=0) return val;
    else return -val;
}

template<typename T>
struct myPoint {
    vector<T> data;

    explicit myPoint(const vector<T> &d){
        data=d;
    }

    myPoint(const myPoint<T> &d){
        data=d.data;
    }

    myPoint<T>& operator=(const myPoint<T> &d){
        data=d.data;
        return *this;
    }

    friend bool operator==(const myPoint<T> &p1, const myPoint<T> &p2){
        size_t size=p1.size();
        double bias=0;
        for(int i=0; i<size;++i){
            bias=max(bias, abs(p1[i]-p2[i]));
        }
        return bias<1e-4;
    }
    T& operator[](std::size_t i){
        return data[i];
    }

    const T& operator[](std::size_t i)const{
        return data[i];
    }

    std::size_t size()const{
        return data.size();
    }
};


struct myPointDoubleHash {
    std::size_t operator()(const myPoint<double> &vec) const {
        std::size_t hash = 14695981039346656037ull;
        for (int i=0; i<vec.size(); ++i) {
            hash ^= (size_t)round(vec[i]*1e4);
            hash *= 1099511628211ull;
        }
        return hash;
    }
};


template<typename T>
struct vectorHash {
    std::size_t operator()(const vector<T> &vec) const {
        std::size_t hash = 14695981039346656037ull;
        for (int i=0; i<vec.size(); ++i) {
            hash ^= (size_t)round(vec[i]*1e4);
            hash *= 1099511628211ull;
        }
        return hash;
    }
};


#endif //UTK_BV_UTK_VECTOR_H
