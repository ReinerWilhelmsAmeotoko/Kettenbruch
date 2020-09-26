//
// Created by reiner on 17/09/2020.
//
// Note to self: To use template class all definitions should go in the .h file
// that is, all the stuff that would usually go in a .cpp file.
// In the cpp file with same class name should be the definitions of
// special functions that are based on a specific type assumption.
//
#include<array>
#include "rationalQ.h"

template<typename Ts>
Ts rationalQ<Ts>::LCM(Ts aa, Ts bb) const  {
    Ts a = labs(aa);
    Ts b = labs(bb);
    if (a == b) {
        return a;
    }
    Ts g = GCD(a, b);
   // std::cout << " GCD of " << aa << " and " << bb << " is " << g << std::endl;
    Ts lc;
    if (g == 1) {
        lc = a * b;
    } else {
        Ts am = a / g;
        Ts bm = b / g;
        lc = am * bm * g;
    }
    //std::cout << " LCM of " << aa << " and " << bb << " is " << lc << std::endl;
    return lc;
}


template<typename Ts>
Ts rationalQ<Ts>::GCD(Ts aa, Ts bb) const {
    Ts a = labs(aa);
    Ts b = labs(bb);
    Ts t;
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
    //std::cout << "GCD(" << aa << "," << bb << ")= " << a << std::endl;
    return a;
}


int findBiggestInArray(const std::array<long int, 25> &coefs, int nmax) {
    long int biggest = coefs[0];
    int index = 0;
    for (int j = 0; j < nmax; ++j) {
        if (coefs[j] >= biggest) {
            index = j;
            biggest = coefs[j];
        }
    }
    return index;
}

/*
 * These functions for continued fractions work with a fixed sized array of long int.
 * I don't like this, as I have to hard program the length in this code.
 * But I leave it for now. It's not that important. I guess I would have to
 * use std::vector instead.
 */
/// computes a fraction of two long integers representing an approximation
/// of the root of a long integer, starting from an appropriately chosen
/// start value (that's where the contention is)
/// \param s      (long int) the number to be taken the sqrt of.
/// \param start  (long int) the starting value for the iteration
/// \param maximalIterations   (int) number of iterations to do. Reasonable is values up to 5.
/// \return  A rational number, that is a fraction of two long integers.

rationalQ<long int> squareRootLong(long int s, long int start, int maximalIterations) {
    rationalQ<long int> S(s);
    rationalQ<long int> sn(start);
    rationalQ<long int> half(1, 2);
    for (unsigned int n = 0; n < maximalIterations; ++n) {
        sn = (sn + S / sn) * half;
        std::cout << n << " iteration result " << sn << std::endl;
    }
    return sn;
}

///
/// \param Frac  input: An rational number as fraction using long int.
/// \param qs    output: an array of long int coefficients (needs to be long enough).
/// \return    nc:  length of the continued fraction, i.e., length of output.
int EuclideanChainLong(const rationalQ<long int> &Frac, std::array<long int, 25> &qs) {
    long int a, b, t, q;
    long int A = Frac.num;
    long int B = Frac.den;
    if (B == 1) {
        qs[0] = A;
        return 1;
    }
    int n = 0;
    a = A;
    b = B;
    while (b != 0 && n < 20) {
        t = b;
        q = a / b;
        qs[n] = q;
        b = a - q * b;
        a = t;
        n += 1;
    }
    //std::cout << "Euclid: n= " << n << std::endl;
    return n;
}

/// This gets the coefficients of a continued fraction and starts at the
/// end to contract it to a simple fraction which is returned.
/// \param n   number of coefficients in qs used.
/// \param qs  and array of long integers containing the coefficients
/// \return    a rational number with long int numerator and denominator
rationalQ<long int> reverseReconstructCF(int n, std::array<long int, 25> &qs) {
    rationalQ<long int> CC(0);
    rationalQ<long int> B(0);
    rationalQ<long int> q(0);
    int i = n;
    while (i > 0) {
        i--;
        q = qs[i];
        CC = B + q;
        B = CC.inverse();
        std::cout << "qs[" << i << "]= " << qs[i] << "  CC= " << CC << " f.p: " << CC.tofloat() << std::endl;
    }
    return CC;
}

double contractFracPartial(int nmax, std::array<long int, 25> &qs) {
    double xn, xnp1, yn, ynp1, anp1, xnm1, ynm1;

    xnm1 = (double) qs[0];
    ynm1 = 1.0;
    xn = (double) qs[0] * (double) qs[1] + 1.0;
    yn = (double) qs[1];
    for (int n = 1; n < nmax; ++n) {
        anp1 = (double) qs[n + 1];
        xnp1 = anp1 * xn + xnm1;
        ynp1 = anp1 * yn + ynm1;
        xnm1 = xn;
        xn = xnp1;
        ynm1 = yn;
        yn = ynp1;
        double tt = xn / yn;
        //std::cout << n << " " << tt << std::endl;
    }
    return xn / yn;
}

bool testPrecision(const rationalQ<long int> sn, rationalQ<long int> S, const double precision) {
    auto exact = sn.tofloat();
    exact *= exact;
    auto Sfl = S.tofloat();
    return std::abs(exact - Sfl) < precision;
};

rationalQ<long int> contractFracPartialRational(int nmax, std::array<long int, 25> &qs) {
    long int xn, xnp1, yn, ynp1, anp1, xnm1, ynm1;

    xnm1 = qs[0];
    ynm1 = 1;
    xn = qs[0] * qs[1] + 1;
    yn = qs[1];
    for (int n = 1; n < nmax; ++n) {
        anp1 = qs[n + 1];
        xnp1 = anp1 * xn + xnm1;
        ynp1 = anp1 * yn + ynm1;
        xnm1 = xn;
        xn = xnp1;
        ynm1 = yn;
        yn = ynp1;
        rationalQ<long int> tt(xn, yn);
        //std::cout << "contract n" << n << " coef " << anp1 << " " << tt << " " << tt.tofloat() << std::endl;
    }
    return rationalQ<long int>(xn, yn);
}
///  An experimental function that computes a rational approximation of a square root of some integer
///  If computes for each iteration a continued fraction for each intermediate approximation, and
///  then contracts that chain fraction to a simple fraction, overwriting the previous approximation
///  In all cases so far observed this results in a much simpler rational number
///  And the number of iterations can be about 8
/// \param s      number from which the square root is to be computed long int
/// \param start  start value  long int
/// \param maximalIterations   integer
/// \param qs     std array of long int length 25. (hard coded). contains on return also the CF coefficients.
/// \return   The final fraction approximating the square root of s
///
rationalQ<long int> squareRootComplicated(rationalQ<long int> S, rationalQ<long int> start, int maximalIterations,
                                          std::array<long int, 25> &qs) {
    rationalQ<long int> sn(start);
    rationalQ<long int> fraction(S);
    rationalQ<long int> previous = sn;
    rationalQ<long int> half(1, 2);
    for (unsigned int n = 0; n < maximalIterations; ++n) {
       std::cout << " In goes:    S= " << S << " sn= " << sn << std::endl;
        rationalQ<long int> tmp = S;
        tmp /= sn;
        tmp += sn;
        tmp /= 2;
        sn = tmp;

       // std::cout << " Out comes:  S= " << S << " sn= " << sn << " tmp " << tmp << std::endl;
        sn.simplify();
       std::cout << " Simplify :  S= " << S << " sn= " << sn << std::endl;
        if (testPrecision(sn, S, 5.0e-12)) {
            std::cout << " at Round " << n << " Test: Precision good enough. Schluss und fertig " << std::endl;
            break;
        }
        if (sn.num <= 0) {
            sn = previous;
            fraction = sn;
            std::cout << " Round " << n << " Mr. Arithmetic banged his head. previous sn : " << sn << std::endl;
            break;
        }
        std::cout.precision(16);
        std::cout << " start Euclidean CF with sn=" << sn << " float: " << sn.tofloat() << std::endl;
        int nc = EuclideanChainLong(sn, qs);
        std::cout << "Euclid needed " << nc << " steps CF: [";
        for (int j = 0; j < nc; ++j) std::cout << " " << qs[j];
        std::cout << " ]" << std::endl;
        if (nc > 20) nc = 20;
        int ncstop = findBiggestInArray(qs, nc);
        std::cout << " biggest coef is at " << ncstop << " value " << qs[ncstop] << std::endl;
        if (ncstop == nc - 1) {
            std::cout << " ncstop zu weit " << std::endl;
        }
        if (ncstop > 3) {
            fraction = contractFracPartialRational(ncstop, qs);
            fraction.simplify();
            sn = fraction;
        }
        std::cout << "Round " << n << ": sn " << sn << " Decimal " << sn.tofloat() << " simplified fraction "
                  << fraction << std::endl;
        if (previous == sn) {
            std::cout << " stopped because there is no change " << std::endl;
            break;
        }
        previous = sn;
        double exact = sn.tofloat();
        double Sfloat = S.tofloat();
        exact *= exact;
        if (std::abs(exact - Sfloat) < 1.0e-14) {
            std::cout << " Precision good enough. Schluss und fertig " << std::endl;
            std::cout << " Original: " << Sfloat << " root^2 : " << exact << std::endl;
            break;
        }
    }
    return sn;
}
