
/** @file desolve.h
 *
 *  Interface to the function desolve implemented in desolve.cpp file. */



#ifndef DESOLVE_H_INCLUDED
#define DESOLVE_H_INCLUDED

#include <ginac/ginac.h>

#define F_expansion 4
#define mF_expansion 5
#define FIM 6

#define riccati 7
#define bernouli 8
#define fracbernouli 9
#define general 10
#define odetype1 11
#define odetype2 12
#define odetype3 13
#define odetype4 14
#define odetype5 15
#define odetype6 16
#define odetype7 17
#define jacobiElip024 18
#define jacobiElip123 19


using namespace std;
using namespace GiNaC;

/** travelling wave coordinate **/
extern ex chi, xi;

extern ex NValue;
extern lst twcPhase, paraInDiffSolve, degAcoeff;



int desolve(const ex& diffeq, const lst& dpndt_vars , const int& method, bool test=false);

//This function can check the solutions of the input differential equation.
//Currently, this function is unable to check all solutions reported by GiNaCDE due to some simplification problems.
//I hope that this function can be improved in future.
//Now, I will suggest to use Maple or Mathematica software to check the solutions.
ex checkSolu(const string& diff_equ, const string& solutions, const string& algebraic_solutions="", const string& solutions_conditions="");

/** summing power of each term and the summed powers are collected in n_pow_clt. **/
class find_n_power
{
public:
    unsigned multindic; // It helps to collect each summed power separately.
    exvector n_pow_clt;
    int n_pow(const ex& _expr, const ex& dpndt_var);
    find_n_power(unsigned _multindic):multindic(_multindic){}
    ~find_n_power(){}
} ;

/** finding order of ode **/
class find_ode_order:public map_function
{
public:
    vector<int> order_clt;
    ex operator()(const ex& _e);
    ~find_ode_order(){}
};

int order(const ex& _expr);


/** finding independent variables in diff eq **/
class find_indpndt_var:public map_function
{
public:
    exset indpndt_var;
    ex operator()(const ex& _e);
    ~find_indpndt_var(){}
};


/** traveling wave form **/
class twf:public map_function
{
    std::map<ex, ex, ex_is_less> indpndt_vars_wt_indx;
public:
    twf(std::map<ex, ex, ex_is_less> _indpndt_vars_wt_indx):indpndt_vars_wt_indx(_indpndt_vars_wt_indx){}
    ex operator()(const ex& _e);
    ~twf(){}
};

/** checking type of diff equations **/
int odeType_check(const ex& _ode, const ex& _dpndtVar);

/** returning solution for Bernouli, Riccati diff equations **/
lst firstOrderDiff_solu(const ex& _indpndt_var, const int& _odetype);

/**checking independent vars in diff. equ. **/
class diffIndpndtPrsntTest:map_function
{
public:
    bool isIndpndtPrsnt;
    exset indpndtVarLst;
    diffIndpndtPrsntTest(bool _isIndpndtPrsnt, exset _indpndtVarLst):isIndpndtPrsnt(_isIndpndtPrsnt),indpndtVarLst(_indpndtVarLst){}
    ex operator()(const ex& _e);
};



#endif // DESOLVE_H_INCLUDED
