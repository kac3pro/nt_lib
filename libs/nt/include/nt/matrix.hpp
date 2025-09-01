#include <optional>
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
    /**
     * return the underlying vector
     */
    std::vector<T> inline vec() const noexcept { return data; }

    ModMatrix(T m, std::size_t width, std::size_t height, std::vector<T> data) : m(m), width(width), height(height), data(data)
    {
        if (width * height != data.size())
            throw std::invalid_argument("data size != width*height");
        for (auto &x : this->data)
            x = posMod(x, m);
    }
    ModMatrix(T m, const ModMatrix<T> &a) : m(m), width(a.width), height(a.height), data(a.vec()) {}
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
            size_t pivotRowIndex;
            bool pivotFound = false;
            for (size_t i = j; i < height && !pivotFound; ++i)
            {
                if ((*this)[i][j] != 0)
                {
                    pivotRowIndex = i; // find a row with a nonzero entry in that column
                    pivotFound = true;
                }
            }
            if (!pivotFound)
                continue; // if you can't, continue
            mulRow(invMod((*this)[pivotRowIndex][j], m), pivotRowIndex);
            for (size_t i = 0; i < height; ++i)
            { // for every other row
                if (i == pivotRowIndex)
                    continue;
                addNTimesIthRowtoJthRow(-(*this)[i][j], pivotRowIndex, i); // subtract accordingly
            }
            if (j < height)
                swapRows(j, pivotRowIndex);
        }
        return (*this);
    }
    ModMatrix<T> &appendColumn(std::vector<T> &col)
    {
        if (col.size() != height)
            throw std::invalid_argument("column size mismatch");
        std::vector<T> newData;
        newData.resize((width + 1) * height);
        for (auto &c : col)
            c = posMod(c, m);
        for (size_t i = 0; i < height; ++i)
        {
            std::copy_n(data.begin() + i * width, width, newData.begin() + i * (width + 1));
            newData[i * (width + 1) + width] = col[i];
        }
        data.swap(newData);
        ++width;
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

/**
 * solve a system of linear equations given by a.x=y if there are infinitely many solutions outputs one of them.
 */
template <typename T>
std::optional<std::vector<T>> solveSystem(ModMatrix<T> a, std::vector<T> y)
{
    if (y.size() != a.height)
        throw std::invalid_argument("y size doesn't match");
    if (!isPrime(a.m))
    {
        std::vector<std::pair<long long, long long>> factors = factorMul(a.m);
        std::vector<std::vector<std::pair<T, T>>> solsModPrimePowers(a.width); //[i][j] stores solutions mod factors[j] of ith variable
        for (auto [p, k] : factors)
        {
            ModMatrix<T> ap(p, a);
            auto sol = solveSystem(ap, y);
            if (!sol)
                return std::nullopt;
            auto solVal = sol.value();
            for (size_t i = 0; i < solVal.size(); ++i)
            {
                solVal[i] = henselLift(static_cast<long long>(solVal[i]), p, k);
                solsModPrimePowers[i].emplace_back(solVal[i], p);
            }
        }
        std::vector<T> globalSols;
        for (auto &solsXi : solsModPrimePowers)
            globalSols.push_back(CRT(solsXi).first);
        return globalSols;
    }
    a.appendColumn(y);
    a.rowReduce();
    std::vector<T> x(a.width - 1, 0);
    for (size_t i = a.height; i-- > 0;)
    {
        bool foundPivot = false;
        size_t pivot;
        for (size_t j = 0; j < a.width - 1; ++j)
        {
            if (a[i][j] == 0)
                continue;
            if (!foundPivot)
            {
                foundPivot = true;
                pivot = j;
            }
        }
        if (foundPivot)
            x[pivot] = addMod(x[pivot], a[i][a.width - 1], a.m);
        else if (y[i] != 0)
            return std::nullopt;
    }
    return x;
}
template <typename T>
T dot(std::span<const T> x, std::span<const T> y, T m)
{
    if (x.size() != y.size())
        throw std::invalid_argument("dot: sizes don't match");
    T res = 0;
    for (int i = 0; i < x.size(); ++i)
        res = addMod(res, mulMod(x[i], y[i], m), m);
    return res;
}
// convenience overloads for std::vector<T> that forward to the span version
template <typename T>
T dot(const std::vector<T> &x, const std::vector<T> &y)
{
    return dot(std::span<const T>(x), std::span<const T>(y));
}

template <typename T>
T dot(std::span<const T> x, const std::vector<T> &y)
{
    return dot(x, std::span<const T>(y));
}

template <typename T>
T dot(const std::vector<T> &x, std::span<const T> y)
{
    return dot(std::span<const T>(x), y);
}
template <typename T>
std::vector<T> dot(const ModMatrix<T> &a, std::span<const T> x)
{
    if (a.width != x.size())
        throw std::invalid_argument("dot sizes don't match");
    std::vector<T> res;
    for (int i = 0; i < a.height; ++i)
    {
        res.push_back(dot(a[i], x, a.m));
    }
    return res;
}
// overload for std::vector argument
template <typename T>
std::vector<T> dot(const ModMatrix<T> &a, const std::vector<T> &x)
{
    return dot(a, std::span<const T>(x));
}