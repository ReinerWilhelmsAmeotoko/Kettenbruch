#include <iostream>
#include <array>
#include <vector>
#include "RatQ.h"
#include "rationalQ.h"
#include "QuaternionQ.h"
#include <math.h>
#include "Kettenbruch.h"
#include <boost/multiprecision/cpp_int.hpp>

rationalQ<long int> squareRootComplicated(rationalQ<long int> S, rationalQ<long int> start, int maximalIterations,
                                          std::array<long int, 25> &qs);

long int myGCD(long int aa, long int bb) {
    long int a = abs(aa);
    long int b = abs(bb);
    long int t;
    if (a < b) {
        t = a;
        a = b;   // store the larger one
        b = t;   // store the smaller one den
    }
    while (b != 0l) {
        t = b;
        b = a % b;
        a = t;   // a is the gcd when the while loop finishes.
    }
    std::cout << "myGCD(" << aa << "," << bb << ")= " << a << std::endl;
    return a;
}

bool checkBoss() {
    // The start array is (5 + x)/(3+y);
    // whereby x is the CF of sqrt(2)= [1;2 2 2 2 2 ..]
    // and y is the CF of sqrt(7) = [2; cycle (1 1 1 4) ]
    // The result is should be a CF of [1;7 2 1 7 1 1 2 1 7 2 3 1 ...]
    // Test not complete \todo complete the checkBoss test
    std::array<std::array<long, 4>, 2> some ({5l,1l, 0l, 0l ,  3l, 0l, 1l, 0l});
    CF::CFMaker<long> Boss(some);
    Boss.printMyTable();
    Boss.receiveUpdate(CF::fromX, 1);
    //Boss.printMyTable();
    Boss.receiveUpdate(CF::fromY, 2);
    //Boss.printMyTable();
    Boss.receiveUpdate(CF::fromX, 2);
    //Boss.printMyTable();
    Boss.receiveUpdate(CF::fromY, 1);
    //Boss.printMyTable();
    Boss.receiveUpdate(CF::fromZ, 2);
    Boss.printMyTable();
    std::array<std::array<long, 4>, 2> should {5l, 10l,  6l, 12l, -4l, -7l, -6l, -11l};
    std::array<std::array<long, 4>, 2> const &is = Boss.lookupTable();
    return (should == is) ;
    // 6   13    6   13
    // 5   10    6   12
    // 5   10    6   12
    //-4   -7   -6  -11
}
int main() {

    CF::processingContFrac();
}

/*

(    1    ,   2    ,    2    ,    3    )
(    2    ,    10    ,    11    ,    15    )
(    4    ,    13    ,    16    ,    21    )
(    2    ,    10    ,    25    ,    27    )
(    12    ,    15    ,    16    ,    25    )
(    12    ,    16    ,    21    ,    29    )
(    2    ,    5    ,    14    ,    15    )
(    4    ,    8    ,    19    ,    21    )
(    2    ,    7    ,    26    ,    27    )
*/
