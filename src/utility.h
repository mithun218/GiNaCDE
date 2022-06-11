
/** @file utility.h
 *
 *  Interface to some GiNaCDE's usefull utilities implemented in utility.cpp file. */



#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include <ginac/ginac.h>
#include <chrono>
#include "inifcns.h"

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif


//#define factor_all 18
#define Gtolerance GiNaC::pow(10,-10)
//#define GiNaCDE_gui

using namespace std;
using namespace GiNaC;

typedef std::map<ex, exvector, ex_is_less> exmapexvec;
typedef std::set<lst,ex_is_less> exsetlst;

/** symbols used in other files **/
extern const symbol symb_,factSymb_, _FAIL;

/** declaring some numbers **/
extern const ex _ex0, _ex1, _ex2, _ex_1, _ex_2, _ex1_2, _ex_1_2, _ex1_4, _ex_1_4;

/** it measures the beginning time of evaluations **/
extern std::chrono::time_point<std::chrono::system_clock> beginTime;

#ifdef GiNaCDE_gui
#include <gtk/gtk.h>
extern GtkWidget *window, *status_bar;
#endif // GiNaCDE_gui

extern string CurrentPath, filename;

/** Dependency variable collector **/
extern exmapexvec var_depend;

/** The variables store solutions and constraints of NLPDEs **/
extern std::vector<lst> solutionClt; extern lst constraints;

bool has_only_digits(const string s);

vector<string> split (const string &s, char delim);

/** Dependency creator **/
class dependc
{
public:
    void clear(const ex& sym_, const ex& _dpnd)
    {
        if(!var_depend[sym_].empty())
            var_depend[sym_].erase(remove(var_depend[sym_].begin(), var_depend[sym_].end(), _dpnd), var_depend[sym_].end());
    }

    void clear(const ex& sym_)
    {
        if(var_depend.find(sym_) != var_depend.end())
            var_depend.erase(sym_);
    }

    exvector get(const ex& sym_)
    {
        return var_depend[sym_];
    }
    int operator()(const ex&, const lst& );
};

extern dependc depend;


/** Symbol finder **/
class symbol_finderc:public map_function
{
public:
    exset symbols;
    void clear()
    {
        symbols.clear();
    }
    ex operator()(const ex& _expr);
    ~symbol_finderc(){}
};

 extern symbol_finderc symbol_finder;
 exset symbols(const ex&);

/** converting number into rational**/
class dorat:public map_function
{
public:
    bool israt;
    dorat(){}
    void set(){israt = true;}
    ex operator()(const ex& _e);
    ~dorat(){}
};

 /** collecting power of each base from pow argument, excludes similar power**/
class basepow_clt:public map_function
{
private:
    exset exprchk;
public:
    std::map<ex, lst, ex_is_less> basepow;
    const ex _var;
    void clear()
    {
        basepow.clear();
        exprchk.clear();
    }
    basepow_clt(const ex _var_):_var(_var_){}
    ex operator()(const ex& _e);
    ~basepow_clt(){}
};

/** Checking presence of functions with given variable in eqsn **/
class funcpresent:public map_function
{
public:
    bool funcpresence = false,
         varInPow = false;
    exset func, funcwtvar;
    ex _var;
    funcpresent(ex _var_):_var(_var_){}
    ex operator()(const ex& _e);
    ~funcpresent(){}
};

/** Calculating Factor of irrational function without expand. **/
ex Factor(const ex&);

/** calculating gcd of list of expressions **/
ex Gcd(lst _exp);

/** calculating gcd of list of expressions **/
ex Lcm(lst _exp);

  /** replacing I by _symb **/
class replaceI:public map_function
{
public:
    ex operator()(const ex& _e);
    ~replaceI(){}
};

  /** Checking polynomial type in fim **/
class polycheckFim:public map_function
{
public:
    bool polytype;
    ex _var;
    polycheckFim(bool _polytype, ex _var_):polytype(_polytype),_var(_var_){}
    ex operator()(const ex& _e);
    ~polycheckFim(){}
};

/// Collecting powers of a variable in fim///
class powClt:public map_function
{
public:
    exvector powers;
    ex _var;
    powClt(ex _var_):_var(_var_){}
    ex operator()(const ex& _e);
    ~powClt(){}
};

 /** Checking polynomial class**/
class polycheck:public map_function
{
public:
    bool polytype;
    ex _var;
    polycheck(bool _polytype, ex _var_):polytype(_polytype),_var(_var_){}
    ex operator()(const ex& _e);
    ~polycheck(){}
};

/** Checking polynomial function **/

inline bool is_poly(const ex& _expr, const ex& _var)
{
    polycheck polycheck(true,_var);
    polycheck(_expr);

    return polycheck.polytype;
}


/// doing conjugate free ///
class conjuFree:public map_function
{
public:
    ex operator()(const ex& _e);
    ~conjuFree(){}
};

extern conjuFree conjuFreee;

/** collecting coefficients of variables in base of non-integer power. **/

inline ex Collect(const ex& _expr, const ex& _var) // This function has been used in odeType_check function
{                                                  // in desolve.cpp file.
    ex temexpr = _expr;

    const ex C1_=reader("C1_");

    if(is_a<power>(_expr) && _expr.op(0).has(_var) && _expr.op(0).is_polynomial(_var)&& _expr.op(1)==_ex1_2)
        temexpr = pow(collect(C1_*expand(_expr.op(0)), _var),_expr.op(1));
    else if(_expr.is_polynomial(_var))
        temexpr = collect(C1_*expand(_expr), _var);

    return temexpr;
}

/** replacing power having "add" base with generated symbols **/ // used in Factor
class powBaseSubs:public map_function
{
    unsigned j;
    ex expr;
    string str;
public:
    bool isNu;
    size_t addNum;
    exmap exprToSymMap;
    powBaseSubs( unsigned j_): j(j_){exprToSymMap.clear();addNum=0;isNu=false;}
    ex operator()(const ex& _e);
    ~powBaseSubs(){}
};

/** this function substitute generated symbols from exmap **/
ex genSymbSubs(const ex& _e, const exmap& highDegSubsClt);

/** Getting lst of coefficients from all terms where present _var.
isCltPowZero = true allow to get coefficients of _var^0. **/
lst collectAllCoeff(const ex& _expr, const lst& _var, const bool& isCltPowZero, exmap& _varsWtcoeff);

/** Getting numerator and denominator.
 *  This functution determines numer/denom
 *  accurately having fractional power.
 *  Such as: Numer_Denom(1/(x/y)^(1/2)) returns {1,(x/y)^(1/2)}, but, cuurently, Ginac,s
 *  numer_denom does not give correct results.
 *  To avoid wrong results, this function replace all the base with fractional power by generated symbol, then it
 *  determines numer/denom.
 *  **/
ex Numer_Denom(const ex& _expr);



#endif // UTILITY_H_INCLUDED
