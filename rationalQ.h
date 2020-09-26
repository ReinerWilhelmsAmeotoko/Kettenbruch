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

    bool operator==(const rationalQ &B) const {
        return (num == B.num && den == B.den);
    }

    bool operator<(const rationalQ &B) const {
        Ts lc = LCM(den, B.den);
        Ts a = (lc / den) * num;
        Ts b = (lc / B.den) * B.num;
        //std::cout << " operator lc, a < b " << lc << " " << a <<" " << b << std::endl;
        return a < b;
    }

    bool operator<=(const rationalQ &B) const {
        Ts lc = LCM(den, B.den);
        Ts a = (lc / den) * num;
        Ts b = (lc / B.den) * B.num;
        //std::cout << " operator lc, a < b " << lc << " " << a <<" " << b << std::endl;
        return a <= b;
    }

    bool operator>(const rationalQ &B) const {
        Ts lc = LCM(den, B.den);
        Ts a = (lc / den) * num;
        Ts b = (lc / B.den) * B.num;
        //std::cout << " operator lc, a < b " << lc << " " << a <<" " << b << std::endl;
        return a > b;
    }

    bool operator>=(const rationalQ &B) const {

        Ts lc = LCM(den, B.den);
        Ts a = (lc / den) * num;
        Ts b = (lc / B.den) * B.num;
        //std::cout << " operator lc, a < b " << lc << " " << a <<" " << b << std::endl;
        return a >= b;
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
    if (den == 1) return;
    if (num == 0) {
        den = 1;
        return;
    }
    if (num == den) {
        num = 1;
        den = 1;
        return;
    }
    Ts g = GCD(num, den);
    num /= g;
    den /= g;
    if (den < 0) {
        num = -num;
        den = -den;
    }
}

template<typename Ts>
std::ostream &operator<<(std::ostream &os, const rationalQ<Ts> &R) {
    os << R.num;
    if (R.den != 1) os << '/' << R.den;
    return os;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator+=(const rationalQ<Ts> &y) {
    Ts ll = LCM(den, y.den);
    Ts na = (ll / den) * num;
    Ts nb = (ll / y.den) * y.num;
    num = na + nb;
    den = ll;
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator-=(const rationalQ<Ts> &y) {
    Ts ll = rationalQ<Ts>::LCM(den, y.den);
    Ts na = (ll / den) * num;
    Ts nb = (ll / y.den) * y.num;
    num = na - nb;
    den = ll;
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator*=(const rationalQ<Ts> &B) {
    num *= B.num;
    den *= B.den;
    reduce();
    return *this;
}

template<typename Ts>
rationalQ<Ts> &rationalQ<Ts>::operator/=(const rationalQ<Ts> &B) {
    num *= B.den;
    den *= B.num;
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
    result += bb;
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
    rationalQ<Ts> bb(b);
    rationalQ<Ts> result(x.num, x.den);
    result += bb;
    return result;
}

template<typename Ts>
rationalQ<Ts> operator-(Ts b, const rationalQ<Ts> &x) {
    rationalQ<Ts> bb(b);
    rationalQ<Ts> result(x.num, x.den);
    result -= bb;
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
rationalQ<Ts> operator*(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    return rationalQ<Ts>(x.num * y.num, x.den * y.den);
}

template<typename Ts>
rationalQ<Ts> operator/(const rationalQ<Ts> &x, const rationalQ<Ts> &y) {
    return rationalQ<Ts>(x.num * y.den, x.den * y.num);
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
