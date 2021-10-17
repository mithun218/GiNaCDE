
/** @file solve.h
 *
 *  Interface to solve function implemented in solve.cpp. */



#ifndef SOLVE_H_INCLUDED
#define SOLVE_H_INCLUDED

#include <ginac/ginac.h>
#include "utility.h"

using namespace std;
using namespace GiNaC;

/** counting total number of add in an equation**/
class totalAddInEq:public map_function
{

public:
    size_t addNum;
    totalAddInEq(){}
    ex operator()(const ex& _e);
    ~totalAddInEq(){}
};


class solvec
{
    lst soluClt;
    lst one_eq_solu;
    exsetlst solu;
    exset SysEquCltEnv;
    totalAddInEq totalAddInEqV;

    lst quadratic(const ex & equ_, const ex& var_);
    ex sqroot(const ex& _exp);
    int one_eq_solutions(const ex& _equ, const ex& _var);
    ex sysequ_red(const exset& sysequ_, const exset& var_);
    bool isVarPrsnt(const ex& _expr, const exset& _var);
    exsetlst solu_subs(set<lst, ex_is_less> solu_);
    lst varOrderByDeg(const lst& low_var, map<unsigned, ex>& eqDivider, unsigned& eqDividerSz);
    lst polySoluWtAutoDegSelct(const ex& _equ, const ex& _var);

public:
    lst cubics(const ex & equ_, const ex& var_);
    lst Ncubics(const ex & equ_, const ex& var_);
    lst Nquartic(const ex & equ_, const ex& var_);
    exsetlst operator()(const lst & equ_, const lst& var_);
    ~solvec(){}
};

/** solve functions **/

inline exsetlst solve(lst equs_, lst vars_)
{
    solvec solvef;

    return solvef(equs_, vars_);
}

/** replacing the "pow" terms with _var with created symbols **/
class powBaseSubsWtVar:public map_function
{
    unsigned j;
    ex expr, _var;
    string str;
public:
    exmap exprToSymMapWtVar;
    powBaseSubsWtVar( unsigned j_,ex _var_): j(j_),_var(_var_){exprToSymMapWtVar.clear();}
    ex operator()(const ex& _e);
    ~powBaseSubsWtVar(){}
};

/** replacing the "pow" terms without _var with created symbols **/
class powBaseSubsWtoutVar:public map_function
{
    unsigned j;
    ex expr, _var;
    string str;
public:
    exmap exprToSymMapWtoutVar;
    powBaseSubsWtoutVar( unsigned j_,ex _var_): j(j_),_var(_var_){exprToSymMapWtoutVar.clear();}
    ex operator()(const ex& _e);
    ~powBaseSubsWtoutVar(){}
};

#endif // SOLVE_H_INCLUDED
