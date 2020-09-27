//
// Created by reiner on 17/09/2020.
//

#ifndef RATIONALS_RATIONALQ_H
#define RATIONALS_RATIONALQ_H

#include <iostream>
#include <array>


template<typename Ts>
class rationalQ {
public:
    Ts num, den;

    rationalQ<Ts>() {
        num = 0;
        den = 1;
    }

    explicit rationalQ<Ts>(Ts a) {
        num = a;
        den = 1;
    }

    rationalQ<Ts>(rationalQ<Ts> &B) {
        num = B.num;
        den = B.den;
        reduce();
    }

    rationalQ<Ts>(const rationalQ<Ts> &B) {
        num = B.num;
        den = B.den;
        reduce();
    }

    void simplify() {
        reduce();
    }

    // note: clang says this should be explicit but that's probably wrong.
    rationalQ<Ts>(int x) {
        num = (Ts) x;
        den = 1l;
    }

    rationalQ<Ts>(Ts a, Ts b) {
        num = a;
        den = b;
        reduce();
    }

    rationalQ<Ts> &operator=(Ts val) {
        num = val;
        den = 1;
        return *this;
    }

    rationalQ<Ts> &operator=(const rationalQ<Ts> &B) {
        num = B.num;
        den = B.den;
        reduce();
        return *this;
    }

    rationalQ<Ts> &operator+=(const rationalQ<Ts> &y);

    rationalQ<Ts> &operator-=(const rationalQ<Ts> &y);

    rationalQ<Ts> &operator*=(const rationalQ<Ts> &B);

    rationalQ<Ts> &operator/=(const rationalQ<Ts> &B);

    rationalQ<Ts> operator-();

    rationalQ<Ts> inverse() const;

    rationalQ<Ts> abs() const;

    bool operator==(const rationalQ<Ts> &B) const {
        return (num == B.num && den == B.den);
    }

    bool operator!=(const rationalQ<Ts> &B) const {
        return (num != B.num || den != B.den);
    }

    bool operator<(const rationalQ<Ts> &B) const {
        Ts a = num;
        Ts b = den;
        Ts c = B.num;
        Ts d = B.den;
        Ts g = GCD(b, d);
        if (g > 1l) {
            b /= g;
            d /= g;
        }
        return (a * d) < (c * b);
    }

    bool operator<=(const rationalQ<Ts> &B) const {
        Ts a = num;
        Ts b = den;
        Ts c = B.num;
        Ts d = B.den;
        Ts g = GCD(b, d);
        if (g > 1l) {
            b /= g;
            d /= g;
        }
        return (a * d) <= (c * b);
    }

    bool operator>(const rationalQ<Ts> &B) const {
        Ts a = num;
        Ts b = den;
        Ts c = B.num;
        Ts d = B.den;
        Ts g = GCD(b, d);
        if (g > 1l) {
            b /= g;
            d /= g;
        }
        Ts ad = a * d;    // this can overflow
        Ts bc = b * c;    // this too. and then the rest isn't working
        return (ad > bc);
    }

    bool operator>=(const rationalQ<Ts> &B) const {
        Ts a = num;
        Ts b = den;
        Ts c = B.num;
        Ts d = B.den;
        Ts g = GCD(b, d);
        if (g > 1l) {
            b /= g;
            d /= g;
        }
        return (a * d) >= (b * c);
    }


    double tofloat() const;

protected:
    Ts LCM(Ts aa, Ts bb) const;

    Ts GCD(Ts aa, Ts bb) const;

    void reduce();

};

/*
friend rationalQ<Ts> operator+(const rationalQ<Ts> &x, const rationalQ<Ts> &y);
friend rationalQ<Ts> operator-(const rationalQ<Ts> &x, const rationalQ<Ts> &y);
friend rationalQ<Ts> operator*(const rationalQ<Ts> &x, const rationalQ<Ts> &y);
friend rationalQ<Ts> operator/(const rationalQ<Ts> &x, const rationalQ<Ts> &y);
*/

template<typename Ts>
void rationalQ<Ts>::reduce() {
    if (den == 1l) return;
    if (num == 0l) {
        den = 1l;
        return;
    }
    if (num == den) {
        num = 1l;
        den = 1l;
        return;
    }
    Ts g = GCD(num, den);
    num /= g;
    den /= g;
    if (den < 0l) {
        num = -num;
        den = -den;
    }
}

template<typename Ts>
std::ostream &operator<<(std::ostream &os, const rationalQ<Ts> &R) {
    os << R.num;
    if (R.den != 1l) os << '/' << R.den;
    return os;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator+=(const rationalQ<Ts> &y) {
    Ts a = num;
    Ts b = den;
    Ts c = y.num, d = y.den;
    Ts g = GCD(b, d);
    if (g > 1l) {
        b /= g;
        d /= g;
    }
    num = a * d + c * b;
    den = b * d * g;
    g = GCD(num, den);
    if (g > 1l) {
        num /= g;
        den /= g;
    }
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator-=(const rationalQ<Ts> &y) {
    Ts a = num;
    Ts b = den;
    Ts c = y.num, d = y.den;
    Ts g = GCD(b, d);
    if (g > 1l) {
        b /= g;
        d /= g;
    }
    num = a * d - c * b; // only difference from +
    den = b * d * g;
    g = GCD(num, den);
    if (g > 1l) {
        num /= g;
        den /= g;
    }
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator*=(const rationalQ<Ts> &B) {
    // reduce size by eliminating possible GCD'S
    // before carrying out the multiplication.
    Ts a = num;
    Ts b = den;
    Ts c = B.num, d = B.den;
    Ts g;
    g = GCD(a, d);
    if (g > 1l) {
        a /= g;
        d /= g;
    }
    g = GCD(b, c);
    if (g > 1l) {
        b /= g;
        c /= g;
    }
    num = a * c;
    den = b * d;
//    num *= B.num;
//    den *= B.den;
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator/=(const rationalQ<Ts> &B) {
    Ts a = num;
    Ts b = den;
    Ts c = B.num, d = B.den;
    // we want: a*d / b*c ;
    Ts g;
    g = GCD(a, c);
    if (g > 1l) {
        a /= g;
        c /= g;
    }
    g = GCD(d, b);
    if (g >= 1l) {
        d /= g;
        b /= g;
    }
    num = a * d;
    den = b * c;
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> operator+(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    rationalQ<Ts> result(x.num, x.den);
    result += y;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator-(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    rationalQ<Ts> result(x.num, x.den);
    result -= y;
    return result;
}


template<typename Ts>
rationalQ<Ts> operator*(const rationalQ<Ts> &x, Ts b) {
    rationalQ<Ts> bb(b);
    rationalQ<Ts> result(x.num, x.den);
    result *= bb;
    return result;
}


template<typename Ts>
rationalQ<Ts> operator*(Ts b, const rationalQ<Ts> &x) {
    rationalQ<Ts> bb(b);
    rationalQ<Ts> result(x.num, x.den);
    result *= bb;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator+(Ts b, const rationalQ<Ts> &x) {
    rationalQ<Ts> result(b);
    result += x;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator-(Ts b, const rationalQ<Ts> &x) {
    rationalQ<Ts> result(b);
    result -= x;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator/(const rationalQ<Ts> &x, Ts b) {
    rationalQ<Ts> bb(b);
    rationalQ<Ts> result(x.num, x.den);
    result /= bb;
    return result;
}


template<typename Ts>
rationalQ<Ts> operator/(Ts b, const rationalQ<Ts> &x) {
    rationalQ<Ts> result(b);
    result /= x;
    return result;
}


template<typename Ts>
rationalQ<Ts> operator*(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    rationalQ<Ts> result(x);
    result *= y;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator/(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    rationalQ<Ts> result(x);
    result /= y;
    return result;
}

template<typename Ts>
rationalQ<Ts> rationalQ<Ts>::operator-() {
    return rationalQ<Ts>(-num, den);
}

template<typename Ts>
double rationalQ<Ts>::tofloat() const {
    return (double) num / ((double) den);
}

template<typename Ts>
rationalQ<Ts> rationalQ<Ts>::inverse() const {
    return rationalQ<Ts>(den, num);
}

template<typename Ts>
rationalQ<Ts> rationalQ<Ts>::abs() const {
    if (num >= 0l) return rationalQ<Ts>(num, den);
    else return rationalQ<Ts>(-num, den);
}

rationalQ<long int> squareRootLong(long int s, long int start, int maximalIterations);

int EuclideanChainLong(const rationalQ<long int> &Frac, std::array<long int, 25> &qs);

rationalQ<long int> reverseReconstructCF(int n, std::array<long int, 25> &qs);

double contractFracPartial(int nmax, std::array<long int, 25> &qs);

rationalQ<long int> contractFracPartialRational(int nmax, std::array<long int, 25> &qs);

template<typename Ts>
rationalQ<Ts>
squareRootComplicated(rationalQ<Ts> S, rationalQ<Ts> start, int maximalIterations, std::array<Ts, 25> &qs);

//rationalQ<long int> squareRootComplicated(rationalQ<long int> s, rationalQ<long int> start, int maximalIterations, std::array<long int, 25> &qs);

#endif //RATIONALS_RATIONALQ_H
