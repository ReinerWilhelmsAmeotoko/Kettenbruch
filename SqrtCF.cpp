//
// Created by reiner on 30/09/2020.
//

#include "SqrtCF.h"
namespace CF {
    void computeCFofSquareRoot(const long A, const long B) {
        CF::SqrtCF<long> mysqrt(A, B);

        std::cout << " Computing the CF of the sqrt(" << A <<"/" << B << ")" << std::endl;
        mysqrt.printCurrentCoefficients();

        // First compute 8 coefficients that aren't stored in the object.
        // this is the usual use: with every makeNextState the next
        // coefficient is computed and not stored.
        // now produce 20 more that are stored in the object.
        // Usually one would call this first just to have a fixed
        // These coefficients then are available by means of
        // nextState, which is a function of the base class
        // Kettenbruch.
        // Generate the continuation of the above started CF.
        mysqrt.start();
        int appended  = mysqrt.makeTheCF(100);
        mysqrt.printCurrentCoefficients();
        unsigned int num = mysqrt.numberCoefficients();
        std::cout << " Appended " << appended << " num Coefficients " << num << std::endl;
        std::cout << " stored CF: " << mysqrt << std::endl;
        std::cout << " reading out stored CF: ";
        for (unsigned int k = 0; k < num; ++k) {
            std::cout << k << ":" << mysqrt.nextState() << " ";
        }
        std::cout << std::endl;

        std::cout << " Next: Continue to poll more 50 more coefficients which won't be stored \n";
        for (unsigned int n=0; n<50; ++n) {
            long qi = mysqrt.makeNextState();
            std::cout << qi << " ";
        }
        std::cout << std::endl;
        num = mysqrt.numberCoefficients();
        std::cout << " Unchanged: num stored Coefficients " << num << std::endl;
        std::cout << " stored CF: " << mysqrt << std::endl;
    }
}


/// Since I could not solve the problem with overflow after around 25-30 iterations
/// for CF's that aren't periodic, I decided to use as template typename
/// a thing from boost which provides 128 bit integers.
/// There is only the problem that I can't print such beasts to std::cout
/// So I have to convert to long double just for printing. The following
/// example gets exactly the first 48 CF coefficients correct,
/// and then  no more (the last 2 are 0 but should be 1 2, according to
/// Mathematica, which by some magic can produce as many coefficients
/// as you like (I'm jealous). With 128 bits, this thing peters out
/// at about round 48, when some of the matrix elements have
/// reached values around 10^37, which has log_2 of 122: That's
/// where the rubber hits the road. 
namespace CF {
    void computeCFof4thRoot(const long A, const long B) {
        CF::ForthRootCF<boost::int128_type> forthRoot(A, B);
        const unsigned int nmax = 50;
        //boost::int128_type
        long  coeffs[nmax];
        for (unsigned int n=0;n<nmax;++n) {
            boost::int128_type  pn = forthRoot.makeNextState();
           coeffs[n] = pn;
        }
        for (unsigned int n=0;n<nmax;++n)  std::cout << " " << coeffs[n] ;
        std::cout <<  std::endl;
    }
}