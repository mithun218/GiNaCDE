
/** @file fim.h
 *
 * Interface to the fim.cpp file*/



#ifndef FIM_H_INCLUDED
#define FIM_H_INCLUDED

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

/** degree of a particular variable **/
class degreePartiVar:public map_function
{
    ex partiVar;
public:
    exvector degreeClt;
    degreePartiVar(ex _partiVar):partiVar(_partiVar){}
    ex operator()(const ex& e);
    ~degreePartiVar(){}
};


class fim
{
    /** removing an expression **/
    ex removeex(const ex& _expr1, const ex& _expr2);

public:
    int operator()(const ex diffeq, const ex dpndt_varChng, const ex dpndt_var, const ex indpndt_var, const lst twc,
                   const ex tw_coordi, const lst phase,
                   const ex tw_coordiPhase,  lst variables, stringstream& solutions, const lst& diffDenomSolu, const ex& remainingDiffpart);
    ~fim(){}
};

lst collectAllCoeff(const ex& _expr, const lst& _varconst, const bool& isCltPowZero);

#endif // FIM_H_INCLUDED
