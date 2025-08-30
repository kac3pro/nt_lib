#include "gtest/gtest.h"
#include "nt.hpp"
#include <utility>
const int BILION = 1'000'000'000;
const long long VERY_LARGE = static_cast<long long>(BILION) * BILION / 100;
TEST(BezoutTest, smallPositive)
{
    int a = 12, b = 18;
    auto [s, t] = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 120;
    b = 1800;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 120890;
    b = 18;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 1220340;
    b = 17;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 1299;
    b = 1800;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = BILION;
    b = BILION + 9;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
}

TEST(BezoutTest, smallAny)
{
    int a = -12, b = 18;
    auto [s, t] = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 120;
    b = -1800;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = 120890;
    b = -18;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = -1220340;
    b = -17;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = -1299;
    b = -1800;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
}

TEST(BezoutTest, large)
{
    long long a = VERY_LARGE + 9, b = VERY_LARGE + 11;
    auto [s, t] = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = -VERY_LARGE * 5;
    b = VERY_LARGE;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
    a = VERY_LARGE * 5;
    b = BILION + 9;
    std::tie(s, t) = bezoutCoef(a, b);
    ASSERT_EQ(a * s + b * t, std::gcd(a, b));
}

TEST(GcdTest, list)
{
    std::vector<long long> a = {100, 10, 50, 20, 202020, VERY_LARGE};
    ASSERT_EQ(gcd(a), 10);
    ASSERT_EQ(gcd(PRIMES), 1);
    std::vector<long long> b;
    for (int i = 0; i < 1'000'000; ++i)
        b.push_back(VERY_LARGE);
    ASSERT_EQ(gcd(b), VERY_LARGE);
    for (int i = 0; i < 1'000'000; ++i)
        b[i] = VERY_LARGE * 11 + 11 * i;
    ASSERT_EQ(gcd(b), 11);
}

TEST(LcmTest, 2values)
{
    ASSERT_EQ(lcm(2, 3), 6);
    ASSERT_EQ(lcm(20, 30), 60);
    ASSERT_EQ(lcm(10, BILION), BILION);
    ASSERT_EQ(lcm(BILION + 9, BILION + 9), BILION + 9);
    ASSERT_EQ(lcm(static_cast<__int128_t>(BILION + 11), static_cast<__int128_t>(BILION + 9)),
              static_cast<__int128_t>(BILION + 9) * (BILION + 11));
    ASSERT_ANY_THROW(lcm(BILION + 11, BILION + 9));
}
TEST(LcmTest, list)
{
    std::vector<long long> a = {100, 10, 50, 20, VERY_LARGE};
    ASSERT_EQ(lcm(a), VERY_LARGE);
    ASSERT_ANY_THROW(lcm(PRIMES));
    std::vector<long long> b;
    for (int i = 0; i < 1'000'000; ++i)
        b.push_back(VERY_LARGE);
    ASSERT_EQ(lcm(b), VERY_LARGE);
}

TEST(InvModTest, test) {
    ASSERT_EQ(invMod(1,13), 1);
    ASSERT_ANY_THROW(invMod(2,6));
    ASSERT_EQ(invMod(2,5), 3);
    int inv = invMod(BILION, BILION+9);
    ASSERT_EQ(invMod(inv, BILION+9), BILION);
    ASSERT_EQ(invMod(1LL,VERY_LARGE), 1);
}

TEST(IsPrimeTest, test) {
    ASSERT_TRUE(isPrime(BILION+9));
    ASSERT_TRUE(isPrime(13));
    ASSERT_TRUE(isPrime(13)*BILION);
    ASSERT_TRUE(!isPrime(static_cast<long long>(BILION+9)*(BILION+7)));
}
TEST(FactorTest, test) {
    ASSERT_EQ(factor(120), std::vector<long long>({2,2,2,3,5}));
    ASSERT_EQ(factor(BILION+7), std::vector<long long>({BILION+7}));
    ASSERT_EQ(factor(static_cast<long long>(BILION+7)*(BILION+9)), std::vector<long long>({BILION+7, BILION+9}));
    ASSERT_EQ(factor(BILION+11), std::vector<long long>({3, 29, 11494253}));
}

template <typename T>
void makeCRTTest(const std::vector<std::pair<T,T>> eqs) {
    auto [x, m] = CRT(eqs);
    for (auto [ai, mi] : eqs) {
        ASSERT_EQ(posMod(x, mi), posMod(ai, mi));
    }   
}
TEST(CRTTest, test) {
    std::vector<std::pair<int,int>> eqs = {{1,3}, {2,5}};
    makeCRTTest(eqs);
    makeCRTTest(std::vector<std::pair<long long, long long>>({{1,3}, {100, BILION+7}, {15,17}}));
    makeCRTTest(std::vector<std::pair<long long, long long>>({{1, BILION+7}, {1,BILION+9}}));
}