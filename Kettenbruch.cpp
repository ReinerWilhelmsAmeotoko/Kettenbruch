//
// Created by reiner on 23/09/2020.
//

#include "Kettenbruch.h"
/// An implementation of Bill Gosper's method to extract continued fractions coefficients from
/// a function z(x,y), where x and y are given in CF coefficients.
/// I made good use of the very good explanations by Mark Jason Dominus on the site
/// https://perl.plover.com/yak/cftalk/
/// The code below is a test example for the use of two template classes in Kettenbruch.h.
/// (described there).
/// An example is the following, which corresponds to the data structure named biggest in the code.
///           27+7x+15y-xy
/// z(x,y) = --------------    where x = sqrt(2) and y = sqrt(7)  is in biggest below.
///           15+3x+10y+2xy
/// In this example I use sqrt(2)=[1; 2 2 2 2 2 . repeated 2.]
/// and sqrt(7) = [2; 1 1 1 4  1 1 1 4 .. repeated 1 1 1 4...]
///
/// Another example is to use data structure xtimesy in initializing CFMaker.
/// This results in sqrt(14) = [3 ; 1 2 1 6  1 2 1 6 repeated ]
///
namespace CF {
    void processingContFrac() {
        std::array<std::array<long, 4>, 2> biggest({27l, 7l, 15l, -1l, 15l, 3l, 10l, 2l});
        std::array<std::array<long, 4>, 2> invertbiggest({15l, 3l, 10l, 2l, 27l, 7l, 15l, -1l});
        std::array<std::array<long, 4>, 2> negated({-27l, -7l, -15l, 1l, 15l, 3l, 10l, 2l});
        std::array<std::array<long, 4>, 2> invertednegated({15l, 3l, 10l, 2l, -27l, -7l, -15l, 1l});
        std::array<std::array<long, 4>, 2> easier({27l, -7l, -5l, 1l, 15l, 3l, 10l, 2l});
        std::array<std::array<long, 4>, 2> xplusy({0l, 1l, 1l, 0l, 1l, 0l, 0l, 0l});
        std::array<std::array<long, 4>, 2> xtimesy({0l, 0l, 0l, 1l, 1l, 0l, 0l, 0l});
        std::array<std::array<long, 4>, 2> idiottest({0l, 0l, 0l, -24l, 0l, 0l, 3l, 0l});

        CF::CFMaker<long> Z(biggest);
        std::string formula = Z.makeFormulaFromMatrix();

        // hard coded CF's of some square roots, used as input. \todo replace by root finding CF
        // The length of these limits the number of recursion steps.
        // X = Sqrt(2): (specifying std::vector<long> is not necessary as seen for Y below)
        CF::Kettenbruch<long> X(std::vector<long>(
                {1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                 2, 2}));
        // Y = Sqrt(7):
        CF::Kettenbruch<long> Y(
                {2, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1,
                 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1,
                 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 4, 1,
                 1, 1});

        std::vector<long> justfortype;          // Empty vector - Trick to create an empty Kettenbruch object
        Kettenbruch<long int> ZZ(justfortype);  // newly generated CF coefficients are stored in this object.
        X.start();   // bring in initial state.
        Y.start();
        long newCoefficient = 0l;                   // When Z can emit a CF coefficients and update it's stored here.
        unsigned int xmax = X.numberCoefficients();
        unsigned int ymax = Y.numberCoefficients();
        unsigned int countX = 0, countY = 0, countZ = 0;
        unsigned int countScaling = 0;
        unsigned int countReducing = 0;
        unsigned int countMultiples = 0;
        unsigned int maxMultiples = 0;
        unsigned int firstScaling = 0;
        unsigned int firstScalingZcount = 0;
        std::cout << " X has " << xmax << " and Y has " << ymax << " coefficients\n";

        // First prime this a little so that X and Y are ahead of Z:
        std::string laterformula = Z.makeFormulaFromMatrix();
        std::cout << " Start Formula: " << laterformula << std::endl;
        for (unsigned int nc = 0; nc < 10; ++nc) {
            Z.receiveUpdate(CF::fromX, X.nextState());
            countX += 1;
            Z.receiveUpdate(CF::fromY, Y.nextState());
            countY += 1;
            if (nc<5) {
                laterformula = Z.makeFormulaFromMatrix();
                std::cout << " Step " << nc + 1 << " " << laterformula << " 2x4 matrix: " << std::endl;
                Z.printMyTable();
            }
        }

        if (Z.tryReducing()) {
            std::cout << "Reducing after Z operated " << countZ << std::endl;
            Z.printMyTable();
        }
        std::cout << " Z's Matrix Before Loop \n";
        Z.printMyTable();

        /// using more than 140 rounds usually doesn't work because of limitations in 64 bit integers.
        /// The other limit is simply how many CF coefficients have been specified for
        /// X and Y. Using 200 here should result in the termination of the algorithm
        /// as X or Y are at the end of their lists.
        for (int kstate = 0; kstate < 200; ++kstate) {  //std::min(xmax, ymax);
            if (Z.tryReducing()) {
                std::cout << "Reducing before Z operated in loop " << countZ << std::endl;
                Z.printMyTable();
                countReducing += 1;
            }
            if (countX > (xmax - 1) || countY > (ymax - 1)) {
                std::cout << " Limit of the coefficients of X or Y reached.\n";
                break;
            }
            unsigned int countZEmit = 0;
            for (int kk = 0; kk < 4; ++kk) {
                if (Z.checkIfZcanEmit(newCoefficient)) {
                    Z.receiveUpdate(fromZ, newCoefficient);
                    countZEmit += 1;
                    ZZ.appendCoeff(newCoefficient);
                    countZ += 1;
                    long theCommonFactor = Z.greatestCommonDivisor();
                    if (theCommonFactor > 1l) std::cout << " all GCD " << theCommonFactor << std::endl;
                    if (Z.tryReducing()) {
                        // this tries to find any GCD in the fractions in Z and divides it out.
                        // So far this has never occurred.
                        std::cout << "Reducing after Z operated in inner loop " << countZ << std::endl;
                        Z.printMyTable();
                        countReducing += 1;
                    }
                } else {
                    if (countZEmit > maxMultiples) maxMultiples = countZEmit;
                    if (kk > 0) countMultiples += 1;
                    break;
                }
            }
            auto tag = Z.checkWhoShouldRun();    // returns an enum tag in CF::emitted.
            switch (tag) {
                case none:     // Z can't decide so we call them both
                    Z.receiveUpdate(fromX, X.nextState());
                    Z.receiveUpdate(fromY, Y.nextState());
                    countX += 1;
                    countY += 1;
                    break;
                case fromX:
                    //std::cout << " calling X \n";
                    Z.receiveUpdate(fromX, X.nextState());
                    countX += 1;
                    break;
                case fromY:
                    //std::cout << " calling Y \n";
                    Z.receiveUpdate(fromY, Y.nextState());
                    countY += 1;
                    break;
            }
            if (kstate > 50) {
                if (Z.tryScaling()) {
                    // This tests if the coefficients are bigger than 2^50 or 2^60 at maximum, and then divides by 4.
                    countScaling += 1;
                    if (firstScaling == 0) firstScaling = kstate + 1;
                    if (firstScalingZcount == 0) firstScalingZcount = countZ;
                }
            }
        }  // end of big loop.

        std::cout << "Z's matrix at the end \n";
        Z.printMyTable();
        std::cout << " ---------- Summary -----------\n";
        std::cout << " X has " << xmax << " and Y has " << ymax << " coefficients\n";
        //std::cout << " X was called " << countX << " times and Y was called " << countY << " times called\n";
        std::cout << " X was called " << X.timesCalled() << " times, and Y was called " << Y.timesCalled()
                  << " times\n";
        std::cout << " Z emitted " << countZ << " coefficients.";
        std::cout << " Z could emit multiples " << countMultiples << " times. Maximal loop: " << maxMultiples
                  << std::endl;
        std::cout << " Scaling was done " << countScaling << " times." << std::endl;
        std::cout << " 1st time scaling round " << firstScaling << " after " << firstScalingZcount << " Z calls "
                  << std::endl;
        std::cout << " Reducing was possible " << countReducing << " times" << std::endl;
        std::cout << " For the function " << formula ;
        std::cout << ", Z emitted " << countZ << " CF coefficients stored in ZZ. \n";
        std::cout << "CF in ZZ: " << ZZ << std::endl;
    }
}