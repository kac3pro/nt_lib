#include <vector>
#include "nt.hpp"

template <typename T>
inline T addMod(T a, T b, T m) { return posMod(a + b, m); }

template <typename T>
inline T mulMod(T a, T b, T m) { return posMod(static_cast<__int128_t>(a) * b, static_cast<__int128_t>(m)); }
template <typename T>
inline T divMod(T a, T b, T m) { return posMod(static_cast<__int128_t>(a) * invMod(b, m), static_cast<__int128_t>(m)); }
template <typename T>
class ModMatrix
{
private:
    std::vector<T> data;

public:
    T m;
    std::size_t width, height;

    ModMatrix(T m, std::size_t width, std::size_t height, std::vector<T> data) : m(m), width(width), height(height), data(data)
    {
        if (width * height != data.size())
            throw std::invalid_argument("data size != width*height");
        for (auto &x : this->data)
            x = posMod(x, m);
    }
    std::span<T> operator[](std::size_t index)
    {
        if (index >= height)
            throw std::out_of_range("row");
        return {data.data() + index * width, width};
    }
    std::span<const T> operator[](std::size_t index) const
    {
        if (index >= height)
            throw std::out_of_range("row");
        return {data.data() + index * width, width};
    }
    ModMatrix<T> &swapRows(size_t i, size_t j)
    {
        if (i >= height || j >= height)
            throw std::out_of_range("row index out of range");
        if (i == j)
            return *this; // nothing to do
        auto it1 = data.begin() + i * width;
        auto it2 = data.begin() + j * width;

        std::swap_ranges(it1, it1 + width, it2);
        return *this;
    }
    ModMatrix<T> &addNTimesIthRowtoJthRow(T n, size_t i, size_t j)
    {
        const auto rowi = (*this)[i];
        auto rowj = (*this)[j];
        for (size_t k = 0; k < rowi.size(); ++k)
        {
            auto &ai = rowi[k];
            auto &aj = rowj[k];
            aj = addMod(aj, mulMod(n, ai, m), m);
        }
        return (*this);
    }
    ModMatrix<T> &mulRow(T n, size_t i)
    {
        auto rowi = (*this)[i];
        for (auto &ai : rowi)
            ai = mulMod(ai, n, m);
        return (*this);
    }

    ModMatrix<T> &rowReduce()
    {
        if (!isPrime(m))
            throw std::runtime_error("rowReduce works only for prime moduli");
        for (size_t j = 0; j < width; ++j)
        { // for each column
            size_t pivotRowIndex = -1;
            for (size_t i = j; i < height; ++i)
            {
                if ((*this)[i][j] != 0)
                    pivotRowIndex = i; // find a row with a nonzero entry in that column
            }
            if (pivotRowIndex == -1) continue; // if you can't, continue
            mulRow(invMod((*this)[pivotRowIndex][j], m), pivotRowIndex);
            for (size_t i = 0; i < height; ++i)
            { // for every other row
                if (i == pivotRowIndex)
                    continue;
                addNTimesIthRowtoJthRow(-(*this)[i][j], pivotRowIndex, i);// subtract accordingly
            }
            if (j < height)
                swapRows(j, pivotRowIndex);
        }
        return (*this);
    }
    static ModMatrix<T> I(size_t d, T m)
    {
        ModMatrix<T> a(m, d, d, std::vector<T>(d * d, 0));
        for (size_t i = 0; i < d; ++i)
        {
            a[i][i] = 1;
        }
        return a;
    }
};
template <typename T>
std::ostream &operator<<(std::ostream &os, ModMatrix<T> const &m)
{
    const std::size_t maxShow = 4;
    os << "ModMatrix(" << m.width << "x" << m.height << ") " << std::endl;
    for (size_t i = 0; i < m.height; ++i)
    {
        for (size_t j = 0; j < m.width; ++j)
        {
            os << m[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}
template <typename T>
inline bool isSameSize(const ModMatrix<T> &left, const ModMatrix<T> &right)
{
    return left.height == right.height && left.width == right.width;
}
template <typename T>
inline bool canAdd(const ModMatrix<T> &left, const ModMatrix<T> &right)
{
    return left.m == right.m && isSameSize(left, right);
}
template <typename T>
inline bool canMultiply(const ModMatrix<T> &left, const ModMatrix<T> &right)
{
    return left.m == right.m && left.width == right.height;
}
template <typename T>
ModMatrix<T> operator+(const ModMatrix<T> &left, const ModMatrix<T> &right)
{
    if (!isSameSize(left, right))
        throw std::invalid_argument("sizes don't match");
    size_t size = left.width * left.height;
    ModMatrix<T> a(left.m, left.width, left.height, std::vector<T>(size, 0));
    for (int i = 0; i < left.height; ++i)
        for (size_t j = 0; j < left.width; ++j)
            a[i][j] = left[i][j] + right[i][j];
    return a;
}
template <typename T>
bool operator==(const ModMatrix<T> &left, const ModMatrix<T> &right)
{
    if (!isSameSize(left, right) || left.m != right.m)
        return false;
    for (size_t i = 0; i < left.height; ++i)
    {
        for (size_t j = 0; j < left.width; ++j)
        {
            if ((left[i][j] - right[i][j]) % left.m != 0)
                return false;
        }
    }
    return true;
}