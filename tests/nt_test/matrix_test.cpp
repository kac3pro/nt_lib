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
TEST(Identity, test) {
    ASSERT_EQ(b, ModMatrix<long long>::I(3,5));
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
    ASSERT_EQ(a[0][2], 1);
}
TEST(RowReduce, small) {
    ASSERT_TRUE(verifyRowReduced(a.rowReduce()));
    ASSERT_TRUE(verifyRowReduced(b.rowReduce()));
    ASSERT_TRUE(verifyRowReduced(c.rowReduce()));
    ASSERT_EQ(d.rowReduce(), ModMatrix<int>::I(3, 2137));
}