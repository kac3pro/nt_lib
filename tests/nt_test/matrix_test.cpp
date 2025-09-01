#include "gtest/gtest.h"
#include "nt/matrix.hpp"
TEST(TestAdd, test) {
    ASSERT_EQ(ModMatrix(2, 2, 2, {1,0,0,1}) + ModMatrix(2,2,2,{1,1,1,1}), ModMatrix(2, 2, 2, {0,1,1,0}));
}

ModMatrix a(13,3,3, {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9 
});
ModMatrix a1(5,3,3, {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9 
});
ModMatrix b(1'000'000'007LL,3,2, {
    2, 4, 1,
    4, 0, 1,
});
ModMatrix c(13,2,3, {
    1, 2, 
    4, 5,
    7, 8,
});
ModMatrix d(2137,3,3, {
    1, 2, 3,
    4, 5, 6,
    7, 8, 11 
});

template <typename T>
bool verifyRowReduced(const ModMatrix<T> &a) {
    int previousPivot = -1;
    for (size_t i = 0; i < a.height; ++i) {
        for (int j = 0; j < a.width; ++j) {
            if (a[i][j] != 0) {
                if(a[i][j] != 1  || j <= previousPivot )
                    return false;
                else {
                    previousPivot = j;
                    break;
                }
            }
            // if (j == a.width - 1)
                // previousPivot  = a.width;
        }
    }
    return true;

}
ModMatrix im(5,3,3, {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1 
});

TEST(Identity, test) {
    ASSERT_EQ(im, ModMatrix<int>::I(3,5));
}
TEST(ElementaryOps, squareBracket) {
    ASSERT_EQ(a[0][0], 1);
    ASSERT_EQ(a[0][1], 2);
    ASSERT_EQ(a[0][2], 3);
    ASSERT_EQ(a[1][0], 4);
}
TEST(ElementaryOps, mulRow) {
    ASSERT_EQ(a.mulRow(2, 0)[0][0], 2);
    ASSERT_EQ(a[0][1], 4);
    ASSERT_EQ(a[0][2], 6);
}
TEST(RowReduce, small) {
    ASSERT_TRUE(verifyRowReduced(a.rowReduce()));
    ASSERT_TRUE(verifyRowReduced(b.rowReduce()));
    ASSERT_TRUE(verifyRowReduced(c.rowReduce()));
    ASSERT_EQ(d.rowReduce(), ModMatrix<int>::I(3, 2137));
}
TEST(Dot, matrix) {
    const std::vector v({1,1,1});
    ASSERT_EQ(dot(a, std::span<const int>(v)), std::vector({6, 2, 11}));
}
TEST (Solve, prime) {
    ModMatrix e5(2, 2,2, {1,1,2,3});
    std::vector<int> y5({1,1});
    auto x5 = solveSystem(e5, y5);
    ASSERT_EQ(y5, dot(e5, x5.value()));
    ModMatrix e1(13,2,2,{1,0,0,1});
    std::vector<int> y({3,2});
    ASSERT_EQ(solveSystem(e1,y), y);
    ModMatrix e2(19, 2,2, {1,1,2,3});
    y.assign({1,1});
    ASSERT_EQ(solveSystem(e2,y), std::vector<int>({2,18}));
    ModMatrix e3(19, 2,2, {1,1,1,1});
    y.assign({1,2});
    ASSERT_EQ(solveSystem(e3,y), std::nullopt);
    ModMatrix e4(7, 2,1, {1, 1});
    y.assign({2});
    auto x = solveSystem(e4,y);
    const std::vector<int> ZERO2({0,0});
    ASSERT_EQ(y,dot(e4,x.value()));
   
}
TEST (Solve, 2by1) {
    ModMatrix e4(143, 2,1, {1, 1});
    std::vector<int> y4({2});
    auto x = solveSystem(e4,y4);
    ASSERT_EQ(y4,dot(e4,x.value()));
}
TEST (Solve, composite) {
    ModMatrix e1(6,2,2,{1,0,0,1});
    std::vector<int> y({3,2});
    auto x1 = solveSystem(e1, y);
    ASSERT_EQ(y, dot(e1, x1.value()));
    ModMatrix e2(30, 2,2, {1,1,2,3});
    y.assign({1,1});
    auto x2 = solveSystem(e2, y);
    ASSERT_EQ(y, dot(e2, x2.value()));
    ModMatrix e3(1'527'955, 2,2, {1,1,1,1});
    y.assign({1,2});
    ASSERT_EQ(solveSystem(e3,y), std::nullopt);
}
