
/** @file integrate.h
 *
 *  Interface to the class integratec and function integrate implemented in integrate.cpp. */



#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

class integratec
{
    const ex var, IntDiffargu;
    ex find_const(const ex expr_);
    ex setup(const ex expr_);
    ex func_inte(ex expr_);
    ex do_inte(const ex expr_);
public:
    integratec(ex var_):var(var_){}
    ex operator()(const ex & expr_);

    ~integratec(){}
};

/** integrating functions **/

inline ex integrate(ex expr_, ex var_)
{
    integratec integratef(var_);

    return integratef(expr_);
}

#endif // INTEGRATE_H_INCLUDED
