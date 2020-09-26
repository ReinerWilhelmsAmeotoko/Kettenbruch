//
// Created by reiner on 16/09/2020.
//

#include "rationalQ.h"
#include "QuaternionQ.h"
#include<array>

//QuaternionQ<long int> sqrtIteration(const QuaternionQ<long int> &Square, const QuaternionQ<long int> &Start, const int maxiter);

QuaternionQ<long>
sqrtIteration(const QuaternionQ<long> &Square, QuaternionQ<long> &Start, const int maxiter);

QuaternionQ<long> sqrtIteration(const QuaternionQ<long> &Square, QuaternionQ<long> &Start, const int maxiter) {
    QuaternionQ<long> xn = Start;
    rationalQ<long> half(1,2);
    for (unsigned int n=0; n<maxiter; ++n) {
          xn = xn + Square/xn;
          xn /= 2;
//        QuaternionQ<long> invxn = xn.inverse();
//        std::cout << " inverse:   n= " << n << " inv(xn)= " <<  invxn << std::endl;
//        QuaternionQ<long> yn = Square*invxn;
//        std::cout << " S*xn^-1    n= " << n << " yn= " << yn << std::endl;
//        xn += yn;
        std::cout << "before *1/2 n= " << n << " yn= " << xn << std::endl;
       // xn /= 2;
        std::cout << "End of loop n= " << n << " xn= " << xn << std::endl;
    }
    return xn;
}

QuaternionQ<long> sqrtScalarIteration(const QuaternionQ<long> &Square, const int maxIterations) {
    std::array<long,25> qs;
    rationalQ<long> absq, div, hh;
    rationalQ<long> one(1l,1l);
    rationalQ<long> nn = Square.Norm();
    std::cout << " First  Sqrt --------------------------- " << std::endl;
    absq = squareRootComplicated(nn, one, maxIterations, qs);
    std::cout << " Sqrt of " << nn << " float: " << nn.tofloat() << " is " <<absq << " float: " << absq.tofloat() << std::endl;
    std::cout << " absq is now " << absq << " float: " << absq.tofloat() << std::endl;
    absq += Square.a0;
    std::cout << " absq is now " << absq << " float: " << absq.tofloat() << std::endl;
    nn = 2l*absq;
    std::cout << " absq is now " << absq << " float: " << absq.tofloat() << std::endl;
    std::cout << " 2nd Sqrt --------------------------- " << std::endl;
    div = squareRootComplicated(nn, one, maxIterations, qs);
    std::cout << " Sqrt of " << nn << " float: " << nn.tofloat() << " is " << div << " float: " << div.tofloat() << std::endl;
    QuaternionQ<long> result;
    result.a0 = absq/div;
    result.vx = Square.vx/div;
    result.vy = Square.vy/div;
    result.vz = Square.vz/div;
    return result;
}
