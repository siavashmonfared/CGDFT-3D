#ifndef CGLINEARSOLVER_H
#define CGLINEARSOLVER_H
typedef double Real;
#include <cmath>
#include <stdlib.h>
#include <iostream>

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// Here ITMAX is the maximum allowed number of iterations, while EPS is a small number to
// rectify the special case of converging to exactly zero function value.
#define ITMAX 5000 // 200, TODO rendre parametrable
#define EPS 1.0e-20 //1.0e-20

#define FREEALL free_fvector(xi,1,n);free_fvector(h,1,n);free_fvector(g,1,n);

// Tolerance passed to brent.
#define TOL 1.0e-6

// Here ITMAX_BRENT is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
// a small number that protects against trying to achieve fractional accuracy for a minimum that
// happens to be exactly zero.
#define ITMAX_BRENT 1000 //100
#define CGOLD 0.3819660
#define ZEPS  1.0e-20 //1.0e-20

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static Real maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

// Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
// maximum magnification allowed for a parabolic-fit step.
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20 // 1.0e-20


Real f1dim (Real x);

Real brent (Real ax, Real bx, Real cx, Real (*f)(Real), Real tol, Real *xmin);
void mnbrak (Real *ax, Real *bx, Real *cx, Real *fa, Real *fb, Real *fc, Real (*func)(Real));
void linmin (Real p[], Real xi[], int n, Real *fret,Real (*func)(Real []));
void frprmn(Real p[], size_t n, Real ftol,size_t fmaxIter, int *iter, Real *fret, Real (*func)(Real []), void (*dfunc)(Real [], Real []));

#endif // CGLINEARSOLVER_H
