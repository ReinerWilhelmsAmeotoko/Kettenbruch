# Kettenbruch
A little project to tinker with continued fractions. (Kettenbruch German: Confinued Fraction). 
An implementation of Bill Gosper's method to extract continued fractions coefficients from
a function z(x,y), where x and y are given in CF coefficients.
I made good use of the very good explanations by Mark Jason Dominus on the site
https://perl.plover.com/yak/cftalk/
Bill Gosper (former MIT) invented a method for programming to create a CF from one or
two others. The algorithm is only described by some examples, and he apparently never
wrote a publication about his method. It's a bit like you have to program it yourself
to understand it. And this is precisely what I did. The code may be "self explanatory":
I put some commentary in it that I hope will be useful. Even though I'm not sure I
will understand in a year from now my own code.
In short Gosper's method relies on statemachines. For example, assume as in my 
programmed example, that I have already the CF for sqrt(2) and for sqrt(7).
Calling sqrt(2) X, and sqrt(7) Y, I want to compute the following function:
z(x,y) =  (27+7x+15y-xy)/(15+3x+10y+2xy)    
where x = sqrt(2) and y = sqrt(7). 
And thereby directly generate the CF for Z(X,Y), of course, all in integer arithmetic.
Gosper's method have then one statemachine Z that can generate the CF of Z(x,y),
while requesting state updates from X and Y. For this simple case, where the
input CF's are hardcoded, e.g., sqrt(2) = [1; 2 2 2 2 2 2 ... ad infinity],
the state update for X and Y are done by simply moving one position forward
and emitting the next coefficient, that is, the sequence of 1, 2, 2, ...
for X and [2, 1 1 1 4 1 1 1 4 ...] for Y, which represents sqrt(7). Z creates its output
CF coeffients by operating on a small coefficient matrix, which is initially
the coefficients in that above formula: 
[[a b c d];[e f g h]] as [[27, 7, 15, -1];[15, 3, 10,  2]]
The updates of this table with new inputs form X and Y follow from 
the algebra used in substiting x --> p+1/x applied to above formula,
and similar for y -- > q + 1/x.  Z updates it's state by replacing it self
as z --> 1/(z - r), where r is an integer number representing the
next coefficient in the output CF. Z has to find the coefficent r by 
checking if current coefficient matrix can be reduced by pulling 
out an integer r, just like one would find r = 1, to rewrite the fraction
23/13 as 1+10/12. Except that in this case the fraction is a low rank rational
function. The value r only depends on the coefficients in the numberator and
denominator polynomials. 
The functions of Z is implemented in the C++ template class called 
CFMaker<typename Ts>, whereas X and Y are represented by a very rudimentary 
template class called Kettenbruch<typename Ts>, which currently is just 
a list or vector of integer numbers. The intention is to also replace 
X and Y by actual generators that spit out the next CF coeffient, that is,
X and Y the can do all the things Z can, but also have a mechanism for
storing and retreiving the CF coefficients. 

I have tested these template classes so far only with a type of long int, that is
64 bit integers (may be 128 on other systems).

I might later write a better description with more details in LaTeX, until then
I highly recommend the https://perl.plover.com/yak/cftalk/ by Mark Jason Dominus
made in 2005. His work and especially of Bill Gosper, the original inventor, are 
hereby happily acknowledged. 
Sept 25, 2020   Reiner Wilhelms-Tricarico.
p.s. This project also contains a class for rational numbers rationalQ, using integer arithmetic
for all usual operations. Furthermore there is a not yet fully developed class for 
Quaternions, which also operate in integer arithmetic, that is, objects of the class 
rationalQ<long> This is still under development.
           
