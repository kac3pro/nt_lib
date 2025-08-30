#include <vector>
#include "nt.hpp"
namespace Matrix {
    
template <typename T>
inline T addMod(T a, T b, T m) { return (a + b) % m;}

template <typename T>
inline T mulMod(T a, T b, T m) { return (static_cast<__int128_t>(a) * b) % m;}
template <typename T>
struct ModMatrix {
    T m;
    std::size_t width;
    std::size_t height;
    std::vector<T> data;

    ModMatrix(T m, std::size_t width, std::size_t height, std::vector<T> data) : 
                    m(m), width(width), height(height), data(std::move(data)) {
                        if (width*height != data.size())
                            throw std::invalid_argument("data size != width*height");
    }
    std::span<T> operator[] (std::size_t index) {
        if (index>=height) throw std::out_of_range("row");
        return {data.data()+index*width, width};
    }
    ModMatrix<T> &switchRows(size_t i, size_t j) {
        if (i >= height || j >= height)
            throw std::out_of_range("row index out of range");
        if (i == j) return *this;                 // nothing to do
        auto it1 = data.begin() + i * width;
        auto it2 = data.begin() + j * width;

        std::swap_ranges(it1, it1 + width, it2);   // swap the entire rows in O(cols)
        return *this;
    }
    ModMatrix<T> &addNTimesIthRowtoJthRow(T n, size_t i, size_t j) {
        auto &rowi = (*this)[i];
        auto &rowj = (*this)[j];
        for (size_t k = 0; k < rowi.size(); ++k) {
            auto &ai = rowi[k];
            auto &aj = rowj[k];
            ai = addMod(ai, mulMod(n, aj, m), m);
        }
        return (*this);
    }
    ModMatrix<T> &mulRow(T n, size_t i) {
        auto &rowi = (*this)[i];
        for (auto &ai : rowi)
            ai = mulMod(ai, n, m);
        return (*this);
    }
    
    ModMatrix<T> &rowReduce() {
        if (!isPrime(m))
            throw std::runtime_error("rowReduce works only form prime moduli");
        
        
        
    }
};
template <typename T>
inline bool isSameSize(const ModMatrix<T> &left, const ModMatrix<T> &right) {
    return left.height == right.height && left.width == right.width;
}
template <typename T>
inline bool canAdd(const ModMatrix<T> &left, const ModMatrix<T> &right) {
    return left.m == right.m && isSameSize(left, right);
}
template <typename T>
inline bool canMultiply(const ModMatrix<T> &left, const ModMatrix<T> &right) {
    return left.m == right.m && left.width == right.height;
}
template <typename T>
ModMatrix<T> operator+(const ModMatrix<T> &left, const ModMatrix<T> &right) {
        if(!isSameSize(left,right))
            throw std::invalid_argument("sizes don't match");
        std::vector<T> arr(left.data.size());
        for (size_t i = 0; i < left.data.size(); ++i) {
            arr[i] = addMod(left.data[i], right.data[i], left.m);
        }
        return ModMatrix<T>(left.m, left.width, left.height, std::move(arr));
    }
}