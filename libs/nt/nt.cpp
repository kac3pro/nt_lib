#include "nt.hpp"
int main() {
    std::vector<long long> res = factor(1'000'000'000'000'000'030LL);
    for (auto a : res) 
        std::cout << a << " ";
    std::cout << "\n";
    return 0;
}