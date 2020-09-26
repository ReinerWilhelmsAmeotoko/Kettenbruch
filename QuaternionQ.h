//
// Created by reiner on 16/09/2020.
//

#ifndef RATIONALS_QUATERNIONQ_H
#define RATIONALS_QUATERNIONQ_H
// A quaternion class based on rational numbers.

#include<iostream>
#include<array>

#include "rationalQ.h"

template<typename Ts>
class QuaternionQ {
public:
    rationalQ<Ts> a0, vx, vy, vz;

    // generators -------
    QuaternionQ() {
        a0 = 0;
        vx = 0;
        vy = 0;
        vz = 0;
    }

    QuaternionQ(rationalQ<Ts> a, rationalQ<Ts> x, rationalQ<Ts> y, rationalQ<Ts> z) {
        a0 = a;
        vx = x;
        vy = y;
        vz = z;
    }

    QuaternionQ(rationalQ<Ts> a) {
        a0 = a;
        vx = 0;
        vy = 0;
        vz = 0;
    }

    QuaternionQ(const std::array<Ts,3> &vect) {
        a0 = 0;
        vx = vect[0];
        vy = vect[1];
        vz = vect[2];
    }

    QuaternionQ(const std::array<rationalQ<Ts> , 3> &vect) {
        a0 = 0;
        vx = vect[0];
        vy = vect[1];
        vz = vect[2];
    }

    // should or should not be explicit?
    QuaternionQ<Ts> &operator=(Ts a) {
        a0 = a;
        vx = 0;
        vy = 0;
        vz = 0;
        return *this;
    }

// -------  generators done.

    rationalQ<Ts> Norm() const;

    QuaternionQ<Ts> conjugate() const;

    QuaternionQ<Ts> inverse() const;
    QuaternionQ<Ts> operator*=(const Ts mul) {
        a0 *= mul;
        vx *= mul;
        vx *= mul;
        vz *= mul;
        return *this;
    }
    QuaternionQ<Ts> operator/=(const Ts mul) {
        a0 /= mul;
        vx /= mul;
        vy /= mul;
        vz /= mul;
        return *this;
    }
    QuaternionQ operator*(const QuaternionQ<Ts> &Y) const;

    QuaternionQ operator/(const QuaternionQ<Ts> &Y) const;

    QuaternionQ operator+(const QuaternionQ<Ts> &Y) const;

    QuaternionQ operator-(const QuaternionQ<Ts> &Y) const;

    QuaternionQ operator+=(const QuaternionQ<Ts> &Y) {
        a0 += Y.a0;
        vx += Y.vx;
        vy += Y.vy;
        vz += Y.vz;
    }
};


template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::conjugate() const {
    return {a0, -vx, -vy, -vz};
}

template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::inverse() const {
    rationalQ<Ts> nn = Norm();
    rationalQ<Ts> v1 = vx/nn;
    rationalQ<Ts> v2 = vy/nn;
    rationalQ<Ts> v3 = vz/nn;
    return {a0 / nn, -v1, -v2, -v3};
}

template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::operator*(const QuaternionQ &Y) const {
    return {
            a0 * Y.a0 - vx * Y.vx - vy * Y.vy - vz * Y.vz,
            vx * Y.a0 + a0 * Y.vx - vz * Y.vy + vy * Y.vz,
            vy * Y.a0 + vz * Y.vx + a0 * Y.vy - vx * Y.vz,
            vz * Y.a0 - vy * Y.vx + vx * Y.vy + a0 * Y.vz
    };
}

template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::operator/(const QuaternionQ &Y) const {
    rationalQ<Ts> nn = Y.Norm();
    return {
            (a0 * Y.a0 + vx * Y.vx + vy * Y.vy + vz * Y.vz) / nn,
            (vx * Y.a0 - a0 * Y.vx + vz * Y.vy - vy * Y.vz) / nn,
            (vy * Y.a0 - vz * Y.vx - a0 * Y.vy + vx * Y.vz) / nn,
            (vz * Y.a0 + vy * Y.vx - vx * Y.vy - a0 * Y.vz) / nn
    };
}

template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::operator+(const QuaternionQ<Ts> &Y) const {
    return {a0 + Y.a0, vx + Y.vx, vy + Y.vy, vz + Y.vz};
}

template<typename Ts>
QuaternionQ<Ts> QuaternionQ<Ts>::operator-(const QuaternionQ<Ts> &Y) const {
    return {a0 - Y.a0, vx - Y.vx, vy - Y.vy, vz - Y.vz};
}

template<typename Ts>
rationalQ<Ts> QuaternionQ<Ts>::Norm() const {
    return rationalQ<Ts>(a0 * a0 + vx * vx + vy * vy + vz * vz);
}

template<typename Ts>
std::ostream &operator<<(std::ostream &os, const QuaternionQ<Ts> &Q) {
    os << " [" << Q.a0 << ", " << Q.vx << ", " << Q.vy << ", " << Q.vz << "]  float:\n";
    os << " [" << Q.a0.tofloat() << ", " << Q.vx.tofloat() << ", " << Q.vy.tofloat() << ", " << Q.vz.tofloat() << "] ";
    return os;
}


QuaternionQ<long int> sqrtIteration(const QuaternionQ<long int> &Square, QuaternionQ<long int> &Start, const int maxiter);
rationalQ<long int> squareRootComplicated(rationalQ<long int> s, rationalQ<long int> start, int maximalIterations, std::array<long int, 25> &qs);

QuaternionQ<long> sqrtScalarIteration(const QuaternionQ<long> &Square, const int maxiterations);

#endif //RATIONALS_QUATERNIONQ_H
