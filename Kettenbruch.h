//
// Created by reiner on 23/09/2020.
//

//#ifndef RATIONALS_KETTENBRUCH_H
//#define RATIONALS_KETTENBRUCH_H


/// Template classes for generating continued fractions from two others.
/// (usually long int is used, and thus far it was only tested for long int)
/// This is built on the original ideas by Bill Gosper, and I also used
/// the very helpful explanations by Mark Jason Dominus on the site
/// https://perl.plover.com/yak/cftalk/
/// There are so far just two template classes:  Kettenbruch<typename Ts> represents merely
/// the coefficients of a continued fraction (CF). The other class, CFMaker<Ts> is the
/// actual machine that generates a new CF for a function z(x,y), see below.
/// The central object stored and manipulated in the class object CFMaker (below) is
/// is a 2x4 matrix of coefficients This matrix is:
/// [[ a b c d ]
///  [ e f g h ]]
/// It is stored in the object of class CFMaker.
/// It is a representation of the coefficients in the function
///           a + bx + cy + dxy
/// z(x,y) = -------------------  and x and y are known continued fractions.
///           e + fx + gy + hxy
/// ------------------------------------ Sept 2020 Reiner Wilhelms-Tricarico. ----
///

#include<iostream>
#include<array>
#include<vector>
#include<string>
#include "rationalQ.h"

#pragma once

namespace CF {
    /// Two enum classes used for control
    /// in the code.
    enum OpType {
        plus, minus, times, divide
    };
    enum emitted {
        fromX, fromY, fromZ, fromXandY, none
    };

/// Kettenbruch<typename>
/// This class is a representation of a continued fraction with somehow computed coefficients.
/// It can be used to directly represent a truncated CF of some expression, e.g., sqrt(2) =[1;2 2 2 2 ...]
/// It is also used to store continued fraction coefficients as they are computed in the class
/// CFMaker below .
/// Name: Kettenbruch is German for continued fraction (chain fraction)

    template<typename Ts>
    class Kettenbruch {
    public:

        Kettenbruch(const std::vector<Ts> &chain) : coefficients(chain) {
            lastIndex = chain.size();
            current = 0;
        }

        Ts Coefficient(unsigned int index) const {
            if (index < coefficients.size()) {
                return coefficients.at(index);
            } else {
                return 0l;
            }
        }
        void countUp() {counting += 1;}
        void start() {
            current = 0;
            counting = 0;
        }

        unsigned int numberCoefficients() const { return coefficients.size(); }

        unsigned int timesCalled() const { return counting; }

        Ts nextState();

        void appendCoeff(const Ts bb) {
            coefficients.push_back(bb);
            lastIndex += 1;
        }

    protected:
        int lastIndex;
        unsigned int current;
        unsigned int counting;
        std::vector<Ts> coefficients;
    };

    template<typename Ts>
    std::ostream &operator<<(std::ostream &os, const Kettenbruch<Ts> &ZZ) {
        os << "[";
        unsigned int num = ZZ.numberCoefficients();
        if (num > 0) {
            os << ZZ.Coefficient(0);
        }
        if (num > 1) {
            os << ";";
        }
        for (unsigned int n = 1; n < num; ++n) os << " " << ZZ.Coefficient(n);
        os << "]";
        return os;
    }

    template<typename Ts>
    Ts Kettenbruch<Ts>::nextState() {
        Ts retval;
        if (current < lastIndex) {
            retval = Coefficient(current);
            current += 1;
        } else {
            retval = 0;
            std::cout << " crossed the ABYSS !!! \n";
        }
        counting += 1;
        return retval;
    }


    template<typename Ts>
    class CFMaker {
    public:
        /// Generator for making a continued fraction from a function computed from two other CF's
        /// The generator expects the entire initial 2x4 matrix of coefficients.
        /// as fixed sized array object.
        ///
        CFMaker(const std::array<std::array<Ts, 4>, 2> &table) : mat2x4(table) {}

        /// Generator that expects one of the four operators in the enum Optype table
        /// and creates the 2x4 coefficient matrix accordingly for plus, minus, times, divide
        ///
        CFMaker(CF::OpType op) { makeInitial2x4matrix(op); }

        /// helper function for debugging. Prints the 2x4 matrix.
        void printMyTable() const {
            for (unsigned int i = 0; i < 2; ++i) {
                for (unsigned int j = 0; j < 4; ++j) std::cout << " " << mat2x4[i][j];
                std::cout << std::endl;
            }
        }

        /// This is to call to update Z (this object) either from emitted signal from X or Y or tell Z to update
        /// its state. The operations correspond to the algebra in (1) fromX:  X-> p+1/X, fromY: Y->p+1/Y
        /// fromXandY:   X->p+1/X and Y->q+1/Y.  Finally fromZ: Z updates to Z -> 1/(Z-p).
        /// this function can be called with 2 or three arguments. The third is only needed
        /// for the operation fromXandY.
        /// In all cases the 2x4 table is updated. See code for details
        Ts receiveUpdate(CF::emitted opcode, Ts p, Ts q = 0l);

        /// calls Z to check if it can emit. If it can update, it returns the emission value
        /// in the argument. But this doesn't update, that's why it's const
        bool checkIfZcanEmit(Ts &emitCoefficient) const;

        /// same as checkIfZcanEmit but with some diagnostic stuff
        bool ifZcanEmit(Ts &emitCoefficient) const;

        /// diagnostic helper function
        int checkoutSituation() const;

        ///  Compares the absolute differences of fractions in the coefficients:
        /// If  |b/f - a/e| > | c/g-a/e| then X should be called next, otherwise Y.
        /// This comparison is implemented by using double floatingpoint to
        /// compare two fractions - otherwise this may not work for
        /// very large numbers.
        CF::emitted checkWhoShouldRun() const;

        /// as soon as the matrix coefficients (a,b, .. h) get too big
        /// this will scale them down by dividing by 4, once they exceed limis
        /// hardcoded in the function (If any coefficient is bigger than 2^60 for example)
        bool tryScaling();

        /// try to find a GCD for the coefficients (a, ..., h). If found all coefficients are
        /// divided by this greatest common divisor.
        bool tryReducing();

        /// computes the GCD(a,b,c,d,e,f,g,h)
        Ts greatestCommonDivisor() const;

        /// Helper function so the caller can inspect the otherwise protected 2x4 table
        const std::array<std::array<Ts, 4>, 2> &lookupTable() const {
            return mat2x4;
        }

        std::string makeFormulaFromMatrix() const;

    protected:
        /// matrix holding coeffients (a b c d; e f g h).
        std::array<std::array<Ts, 4>, 2> mat2x4;

        /// called by one of the generators to define the initial matrix according to the
        /// operations x+y, x-y, x*y, x/y. as given by the enum: plus, minus, times, divide
        void makeInitial2x4matrix(CF::OpType op);

        /// set's all elements olf the 2x4 matrix to zero.
        void clearTwofour() {
            for (unsigned int i = 0; i < 2; ++i) {
                for (unsigned int j = 0; j < 4; ++j) mat2x4[i][j] = 0l;
            }
        }

        /// Least common multiplier of two integer numbers (helper function)
        /// may not be used currently.
        Ts LCM(Ts aa, Ts bb) const;

        /// Greatest Common Divisor or two integer numbers
        Ts GCD(Ts aa, Ts bb) const;
    };

    template<typename Ts>
    Ts CFMaker<Ts>::receiveUpdate(CF::emitted opcode, Ts p, Ts q) {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        switch (opcode) {
            case fromX:   // x-> p+1/x changes the coefficient matrix to:
                mat2x4 = {b, a + b * p, d, c + d * p, f, e + f * p, h, g + h * p};
                break;
            case fromY:   // y-> p+1/y changes the coefficient matrix to:
                mat2x4 = {c, d, a + c * p, b + d * p, g, h, e + g * p, f + h * p};
                break;
            case fromZ:   //  z -> 1/(z-p)
                mat2x4 = {e, f, g, h, a - e * p, b - f * p, c - g * p, d - h * p};
                break;
            case fromXandY:  // x-> p+1/x and y-> q+1/y (requires the 3rd argument)
                mat2x4 = {d, c + d * p, b + d * q, a + b * p + c * q + d * p * q,
                          h, g + h * p, f + h * q, e + f * p + g * q + h * p * q};
                break;
        }
        Ts alle = greatestCommonDivisor();
        //std::cout << " Greatest Common Divisor " << alle << std::endl;
        return p;
    }

    template<typename Ts>
    void CFMaker<Ts>::makeInitial2x4matrix(CF::OpType op) {
        clearTwofour();
        switch (op) {
            case plus:
                mat2x4[0][1] = 1l;
                mat2x4[0][2] = 1l;
                mat2x4[1][0] = 1l;
                break;
            case minus:
                mat2x4[0][0] = 1l;
                mat2x4[0][2] = -1l;
                mat2x4[1][0] = 1l;
                break;
            case times:
                mat2x4[1][0] = 1l;
                mat2x4[0][3] = 1l;
                break;
            case divide:
                mat2x4[0][1] = 1l;
                mat2x4[1][2] = 1l;
                break;
            default:
                break;
        }
    }

    template<typename Ts>
    bool CFMaker<Ts>::checkIfZcanEmit(Ts &emitCoefficient) const {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        Ts r1 = 0l, r2 = 0l, r3 = 0l, r4 = 0l;
        if (e != 0l) r1 = a / e;
        if (f != 0l) r2 = b / f;
        if (g != 0l) r3 = c / g;
        if (h != 0l) r4 = d / h;
        if (r1 == r2 && r2 == r3 && r3 == r4) {
            emitCoefficient = r1;
            return true;
        } else {
            emitCoefficient = 0l;
            return false;
        }
    }

    template<typename Ts>
    CF::emitted CFMaker<Ts>::checkWhoShouldRun() const {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];

        if (e > 0l && f > 0l && g > 0l) {
            rationalQ<Ts> r0(a, e);
            rationalQ<Ts> rx(b, f);
            rationalQ<Ts> ry(c, g);
            double r0f = r0.tofloat();
            double rxf = rx.tofloat();
            double ryf = ry.tofloat();
            double d1 = std::abs(rxf - r0f);
            double d2 = std::abs(ryf - r0f);
            if (d1 > d2) {
                return fromX;
            } else {
                return fromY;
            }
        }
        return none;
    }

    template<typename Ts>
    bool CFMaker<Ts>::tryScaling() {
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

    template<typename Ts>
    bool CFMaker<Ts>::tryReducing() {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        Ts myGCD, g2, g3, g4;
        myGCD = greatestCommonDivisor();
        if (myGCD > 1l) {
            mat2x4 = {a / myGCD, b / myGCD, c / myGCD, d / myGCD, e / myGCD, f / myGCD, g / myGCD, h / myGCD};
            return true;
        }
        return false;
    }

    template<typename Ts>
    Ts CFMaker<Ts>::LCM(Ts aa, Ts bb) const {
        Ts a = labs(aa);
        Ts b = labs(bb);
        if (a == b) {
            return a;
        }
        Ts g = GCD(a, b);
        Ts lc;
        if (g == 1) {
            lc = a * b;
        } else {
            Ts am = a / g;
            Ts bm = b / g;
            lc = am * bm * g;
        }
        return lc;
    }

    template<typename Ts>
    Ts CFMaker<Ts>::GCD(Ts aa, Ts bb) const {
        Ts a = labs(aa);
        Ts b = labs(bb);
        Ts t;
        if (a == 0l || b == 0l) return std::max(a, b);
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
        return a;
    }

    template<typename Ts>
    Ts CFMaker<Ts>::greatestCommonDivisor() const {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        Ts g1, g2, g3, g4, g12, g22, ggall;
        g1 = GCD(a, e);
        g2 = GCD(b, f);
        g3 = GCD(c, g);
        g4 = GCD(d, h);
        g12 = GCD(g1, g2);
        g22 = GCD(g3, g4);
        ggall = GCD(g12, g22);
        return ggall;
    }

    template<typename Ts>
    bool CFMaker<Ts>::ifZcanEmit(Ts &emitCoefficient) const {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        Ts r1 = 0l, r2 = 0l, r3 = 0l, r4 = 0l;
        int situation = checkoutSituation();
        if (situation == -1) std::cout << " negative \n";
        if (situation == 0) {
            std::cout << " indifferent \n";
            printMyTable();
        }
        if (situation == 1) std::cout << "positive \n";
        if (e != 0l) r1 = a / e;
        if (f != 0l) r2 = b / f;
        if (g != 0l) r3 = c / g;
        if (h != 0l) r4 = d / h;
        if (r1 == r2 && r2 == r3 && r3 == r4) {
            emitCoefficient = r1;
            return true;
        } else {
            emitCoefficient = 0l;
            return false;
        }
    }

    template<typename Ts>
    int CFMaker<Ts>::checkoutSituation() const {
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
        bool upperbig =
                ((aa >= ee || ee == 0l) && (bb >= ff || ff == 0l) && (cc >= gg || gg == 0l) && (dd >= hh || hh == 0l));
        bool lowerbig =
                ((ee > aa || aa == 0l) && (ff > bb || bb == 0l) && (gg > cc || cc == 0l) && (hh > dd || dd == 0l));

        bool upperbigsans = ((aa >= ee) && (bb >= ff) && (cc >= gg) && (dd >= hh));
        bool lowerbigsans =
                ((ee > aa) && (ff > bb) && (gg > cc) && (hh > dd));
        if (upperbigsans && lowerbigsans) return 999;
        if (upperbig) return +1;
        if (lowerbig) return -1;
        std::cout << " abs upper values: " << aa << " " << bb << " " << cc << " " << dd << std::endl;
        std::cout << " abs lower values: " << ee << " " << ff << " " << gg << " " << hh << std::endl;
        return 0;
    }

    template<typename Ts>
    std::string CFMaker<Ts>::makeFormulaFromMatrix() const {
        Ts a = mat2x4[0][0], b = mat2x4[0][1], c = mat2x4[0][2], d = mat2x4[0][3];
        Ts e = mat2x4[1][0], f = mat2x4[1][1], g = mat2x4[1][2], h = mat2x4[1][3];
        bool previous = false;
        int countNum = 0, countDen = 0;
        for (unsigned int k = 0; k < 4; k++) {
            if (mat2x4[0][k] != 0l) countNum += 1;
            if (mat2x4[1][k] != 0l) countDen += 1;
        }
        std::string formula = "Z(X,Y)= ";
        if (countNum > 1) formula += "(";
        if (a != 0l) {
            formula += std::to_string(a);
            previous = true;
        }
        if (b != 0l) {
            if (b > 0l && previous) formula += "+" + std::to_string(b);
            else formula += std::to_string(b);
            formula += "X";
            previous = true;
        }
        if (c != 0l) {
            if (c > 0l && previous) formula += "+" + std::to_string(c);
            else formula += std::to_string(c);
            formula += "Y";
            previous = true;
        }
        if (d != 0l) {
            if (d > 0l && previous) formula += "+" + std::to_string(d);
            else formula += std::to_string(d);
            formula += "XY";
            previous = true;
        }
        if (countNum > 1) formula += ")";
        formula += "/";
        previous = false;
        if (countDen > 1) formula += "(";
        if (e != 0l) {
            formula += std::to_string(e);
            previous = true;
        }
        if (f != 0l) {
            if (f > 0l && previous) formula += "+" + std::to_string(f);
            else formula += std::to_string(f);
            formula += "X";
            previous = true;
        }
        if (g != 0l) {
            if (g > 0l && previous) formula += "+" + std::to_string(g);
            else formula += std::to_string(g);
            formula += "Y";
            previous = true;
        }
        if (h != 0l && previous) {
            if (h > 0l && previous) formula += "+" + std::to_string(h);
            else formula += std::to_string(h);
            formula += "XY";
        }
        if (countDen > 1) formula += ")";
        formula += " ";
        return formula;
    }

//// other tokens

    void processingContFrac();


}  //end namespace CF
