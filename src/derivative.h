
/** @file derivative.h
 *
 *  Interface to the derivative.cpp file. */


#ifndef DERIVATIVE_H_INCLUDED
#define DERIVATIVE_H_INCLUDED

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

extern lst twcPhase;
 /** substituting in Int and Diff arguments  and mapping**/
class subs_IntDiffargu:public map_function
{
    unsigned i,j;
    ex _expr, _var, _temexpr;
    exset symclt;
    string str;
public:
    exmap  IntDiffargu_map, Conjurgu_map;
    subs_IntDiffargu(unsigned i_, unsigned j_, ex _var_):i(i_), j(j_), _var(_var_){}
    ex operator()(const ex& _e);
    ~subs_IntDiffargu(){}
};

/** Derivative for funtions with dependency **/
ex pdiff(const ex&, const ex&, const ex&);

/** Evaluation of Diff, Integrate **/
class evaluatec:public map_function
{
    ex _tem;
public:
    ex operator()(const ex& _e);
    ~evaluatec(){}
};
extern evaluatec evaluateg;

ex evaluate(const ex&);

ex zero_order_rem(const ex&);
#endif // DERIVATIVE_H_INCLUDED
