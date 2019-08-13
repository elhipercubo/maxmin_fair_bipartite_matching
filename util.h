#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <numeric>
#include <locale.h>
#include <map>
#include <iostream>
#include <queue>
#include <set>
#include <vector>
#include <unistd.h>
#include <boost/rational.hpp>

#define gc() getchar_unlocked()

int isSpaceChar(char c) {
    return c == ' ' || c == '\n' || c == '\r' || c == '\t' ;
}

// Fast input
inline int read_int(int* ret)
{
    char ch;
    int val=0;
    ch=gc();
    while(isSpaceChar(ch))
        ch=gc();
    if (feof(stdin)) return 0;
    while(!isSpaceChar(ch))
    {
        val=(val*10)+(ch-'0');
        ch=gc();
    }
    *ret = val;
    return 1;
}

typedef boost::rational<long long> Rational;
using boost::rational_cast;
namespace boost {
    template <typename IntType>
    constexpr IntType floor(rational<IntType> const& r) {
        return static_cast<IntType>(r.numerator() / r.denominator());
    }
    template <typename IntType>
    constexpr IntType ceil(rational<IntType> const& r) {
        return static_cast<IntType>((r.numerator() + r.denominator() - 1) / r.denominator());
    }
}
