
/** @file F_expns_method.h
 *
 * Interface to the F_expns_method.cpp file. */



#ifndef F_EXPNS_METHD_H_INCLUDED
#define F_EXPNS_METHD_H_INCLUDED
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

extern const symbol _n;


extern bool ASolve;
extern bool positivePart, negativePart;


class F_expans
{

public:
    int operator()(const ex diffeq, const ex dpndt_varChng, const ex dpndt_var, const ex indpndt_var, const lst twc,
                   const ex tw_coordi, const lst phase, const ex tw_coordiPhase, lst variables, stringstream& solutions,
                   const ex& Nvalue, const int method, const lst& diffDenomSolu, const ex& remainingDiffpart);
    ~F_expans(){}
};

#endif // F_EXPNS_METHD_H_INCLUDED
