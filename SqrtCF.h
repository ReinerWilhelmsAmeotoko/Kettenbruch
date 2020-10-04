//
// Created by reiner on 30/09/2020.
//

#ifndef RATIONALS_SQRTCF_H
#define RATIONALS_SQRTCF_H

#include<iostream>
#include<array>
#include<vector>
#include<string>
#include <boost/multiprecision/cpp_int.hpp>
#include "rationalQ.h"
#include "Kettenbruch.h"

#pragma once

namespace CF {
///
/// \tparam Ts Usually long int
/// SqrtCF is a derived class from Kettenbruch.
/// This class can generate the coefficients of square roots of fractions on demand
/// or stored in an array of the base class.
/// It receives in the Generator a fraction A/B. This is used to set up
/// the machinery for computing CF coefficients. Whenever the coefficient
/// is desired, call the method makeNextState(), which returns a value
/// of the \tparam, usually long.
///
    template<typename Ts>
    class SqrtCF : public Kettenbruch<Ts> {
    public:
        /// the generator has two parameters, namely integer value of
        /// a numerator and denominator of a fixed fraction. For example, Called with
        /// (7,1) it will generate the Cf of sqrt(7), or (23,13) of 23/13.
        SqrtCF(const Ts Num, const Ts Den) : Kettenbruch<Ts>(std::vector<Ts>()) {
            std::cout << " SqrtCF generator called \n";
            init(Num, Den);
        }

        SqrtCF() : Kettenbruch<Ts>(std::vector<Ts>()) {
            std::cout << " SqrtCF generator called sanse arguments\n";
        }

        /// computes niterations CF coefficients in one go, and they are stored
        /// in the the next CF coefficient and update the system matrix.
        /// And then return this to the caller.
        int makeTheCF(const int niterations);

        /// helper function for diagnosis
        void printCurrentCoefficients() const {
            using boost::multiprecision::cpp_int;
            //std::numeric_limits<boost::multiprecision::int128_t> bull00 = mat2x2[0][0];
            //std::numeric_limits<boost::multiprecision::int128_t> bull01 = mat2x2[0][1];
            //std::numeric_limits<boost::multiprecision::int128_t> bull10 = mat2x2[1][0];
            //std::numeric_limits<boost::multiprecision::int128_t> bull11 = mat2x2[1][1];
            std::cout << "(" << (long double) mat2x2[0][0] << " Y + " << (long double) mat2x2[0][1] << ")/("
                      <<  (long double)mat2x2[1][0] << " Y + " << (long double) mat2x2[1][1] << ") " << std::endl;
        }

        /// compute the next CF coefficient and update the system matrix.
        /// And then return this to the caller.
        Ts makeNextState();

        void init(const Ts aa, const Ts bb) {
            A = aa;
            B = bb;
            mat2x2 = {0l, A, B, 0l};
        }

    protected:

        Ts iteratedRoot();

        Ts A, B;
        std::array<std::array<Ts, 2>, 2> mat2x2;
    };

    template<typename Ts>
    int SqrtCF<Ts>::makeTheCF(const int niterations) {
        Ts q = 1l;
        unsigned int nn = 0;
        while (nn < niterations) {
            q = makeNextState();
            Kettenbruch<Ts>::appendCoeff(q);
            nn += 1;
        }
        return (int) nn;
    }

    template<typename Ts>
    Ts SqrtCF<Ts>::iteratedRoot() {
        Ts a = mat2x2[0][0], b = mat2x2[0][1];
        Ts c = mat2x2[1][0], d = mat2x2[1][1];
        double af = a, bf = b, cf = c;
        double y = 1.0;
        Ts retval = 1l;
        if (c == 0l && a != 0l) {
            y = -bf / (2 * af);
            retval = (Ts) y;
            return retval;
        } else {
            y = 2 * af / cf + 1.0;
            for (unsigned int k = 0; k < 10; ++k) {
                double den = cf * y - af;
                if (den < 1.0e-30) {
                    std::cout << " Underflow in c*y-a \n";
                    printCurrentCoefficients();
                    return 1l;
                }
                y = ((af * y + bf) / den + y) / 2;
            }
            retval = (Ts) y;
            return retval;
        }
    }

    template<typename Ts>
    Ts SqrtCF<Ts>::makeNextState() {
        Ts q = iteratedRoot();
        Ts a = mat2x2[0][0];
        Ts b = mat2x2[0][1];
        Ts c = mat2x2[1][0];
        Ts d = mat2x2[1][1];
        mat2x2 = {c * q - a, c, b + a * q + a * q - c * q * q, a - c * q};
        return q;
    }
//------------- End of SqrtCF

/// A class that computes the sqrt(sqrt(Num/Den)) and produces its CF.
/// It relies on a similar object (for now protected), called
/// MisterX who sends the CF coefficients of the simple sqrt(Num/Den);

    template<typename Ts>
    class ForthRootCF : public Kettenbruch<Ts> {
    public:
        ForthRootCF(const Ts Num, const Ts Den) : Kettenbruch<Ts>(std::vector<Ts>()) {
            std::cout << " ForthRootCF generator called \n";
            init(Num, Den);  // initializes also MisterX.
        }

        void init(const Ts aa, const Ts bb) {
            A = aa;
            B = bb;
            mat2x4 = {0l, 1l, 0l, 0l, 0l, 0l, 1l, 0l};   // X/Y
            std::cout << " Calling init on MisterX \n";
            MisterX.init(A, B);  // This guy sends the CF coefficients of a sqrt(aa/bb) on demand.
            numXcalled = 0;
        }

        /// compute the next CF coefficient and update the system matrix.
        /// And then return this to the caller.
        Ts makeNextState();

        void printMyMatrix() const;

    protected:

        void updateEmitX(const Ts q);

        void trippleUpdate(const Ts q, const Ts p);

        Ts iterativeFind(const Ts qn);

        bool tryScaling();

        // -------
        Ts A, B;
        std::array<std::array<Ts, 4>, 2> mat2x4;
        int numXcalled;
        SqrtCF<Ts> MisterX;
    };

    template<typename Ts>
    Ts ForthRootCF<Ts>::makeNextState() {
        Ts qn, pn;
        if (numXcalled <= 0) {
            numXcalled = 0;
            for (unsigned int n = 0; n < 5; ++n) {
                qn = MisterX.makeNextState();
                //std::cout << " MisterX " << n << " " << qn << std::endl;
                //printMyMatrix();
                updateEmitX(qn);
                numXcalled += 1;
            }
        }
        //std::cout << " X called " << numXcalled << std::endl;
        qn = MisterX.makeNextState();
        //std::cout << " MisterX called " << numXcalled << " " << qn << std::endl;
        numXcalled += 1;
        pn = iterativeFind(qn);
        trippleUpdate(qn, pn);
       // if (tryScaling()) {
       //     std::cout << " scaling at " << numXcalled << std::endl;
       // }
        return pn;
    }

    template<typename Ts>
    void ForthRootCF<Ts>::updateEmitX(const Ts q) {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        mat2x4 = {b, a + b * q, d, c + d * q,
                  f, e + f * q, h, g + h * q};
    }

    template<typename Ts>
    Ts ForthRootCF<Ts>::iterativeFind(const Ts qn) {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        long double A, B, C, D, G, H;
        A = a + qn * b;
        B = c + qn * d;
        C = e + qn * f;
        D = g + qn * h;
        G = (2.0 * B) / ((double) D);
        H = A / D;
        long double yn = G;
        if (yn < 2.0) yn=2.0;
        // std::max(G, 1.0);
        std::cout << " A=" << A << " B=" << B << " C=" << C<< " D=" << D<< " G=" << G<< " H=" << H << std::endl;
        std::cout << " " << yn;
        for (unsigned int m = 0; m < 25; ++m) {
            yn = (yn * yn + H) / (2.0 * yn - G);
            std::cout << " " << yn;
        }
        std::cout << std::endl;
        Ts retval = (Ts) yn;
        std::cout << "Iteration result: yn= " << yn << " floor: " << (long double) retval << std::endl;
        return retval;
    }

    template<typename Ts>
    void ForthRootCF<Ts>::trippleUpdate(const Ts q, const Ts p) {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        mat2x4 = {d, c + d * q, b + d * p, a + b * q + p * (c + d * q),
                  h, g + h * q, f + h * p, e + f * q + p * (g + h * q)};
        //printMyMatrix();
        a = mat2x4[0][0];
        b = mat2x4[0][1];
        c = mat2x4[0][2];
        d = mat2x4[0][3];
        e = mat2x4[1][0];
        f = mat2x4[1][1];
        g = mat2x4[1][2];
        h = mat2x4[1][3];
        mat2x4 = {e, f, g, h, a - e * p, b - f * p, c - g * p, d - h * p};
        printMyMatrix();
    }

    template<typename Ts>
    void ForthRootCF<Ts>::printMyMatrix() const {
       // Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
       // Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        for (unsigned int k=0; k<2; ++k){
            for (unsigned int j=0 ; j < 4; ++j) std::cout << (long double) mat2x4[k][j] << " ";
            std::cout << std::endl;
        }

       // std::cout << "[" << a << "," << b << "," << c << "," << d << ";\n "
       //           << e << "," << f << "," << g << "," << h << "]\n";
    }

    template<typename Ts>
    bool ForthRootCF<Ts>::tryScaling() {
        const Ts mx = 1152921504606846976l; // 2^60 top
        const Ts mb = 1125899906842624l;    // 2^50 bottom
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        Ts aa = std::abs(a);
        Ts bb = std::abs(b);
        Ts cc = std::abs(c);
        Ts dd = std::abs(d);
        Ts ee = std::abs(e);
        Ts ff = std::abs(f);
        Ts gg = std::abs(g);
        Ts hh = std::abs(h);
        if (aa > mb && bb > mb && cc > mb && dd > mb && ee > mb && ff > mb && gg > mb && gg > mb) {
            if (aa > mx || bb > mx || cc > mx || dd > mx || ee > mx || ff > mx || gg > mx || hh > mx) {
                //mat2x4 = {a / 2l, b / 2l, c / 2l, d / 2l, e / 2l, f / 2l, g / 2l, h / 2l};
                mat2x4 = {a / 4l, b / 4l, c / 4l, d / 4l, e / 4l, f / 4l, g / 4l, h / 4l};
                return true;
            }
        }
        return false;
    }

    void computeCFofSquareRoot(const long A, const long B);

    void computeCFof4thRoot(const long A, const long B);
} // -- namespace CF

#endif //RATIONALS_SQRTCF_H
