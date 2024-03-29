
/** @file simplify.h
 *
 *  Interface to GiNaCDE's simplify function implemented in simplify.cpp. */


#ifndef SIMPLIFY_H_INCLUDED
#define SIMPLIFY_H_INCLUDED

#define AlgSimp 0
#define TrigSimp 1
#define TrigCombine 2
#define logSimp 3
#define JacobiSimp 4
#define AlgSimp2 5
#define HyperSimp 6
#define FuncSimp 7


#include <iostream>
#include <stdexcept>
#include <ginac/ginac.h>
#include "utility.h"

using namespace std;
using namespace GiNaC;

extern size_t expandLevel, addNumFrFactr;// these variables limit the simplification of algebraic expressions
extern long long int largstNumsimp; // this is the maximum number for simplification

class simplifyc
{
    //int rules = AlgSimp;
    exmap AlgSimpRules1,AlgSimpRules2, AlgSimpRules3, TrigSimpRules1, TrigSimpRules2, HyperSimpRules1, HyperSimpRules2,
          TrigCombineRules, logSimpRules, JacobiSimpRules1, JacobiSimpRules2;
    int SetRules(int m = AlgSimp);
public:
    simplifyc(){}
    ex operator()(const ex& e,const int& rules = AlgSimp);
    ~simplifyc(){}
};

/**This collect all common factors including numerical numbers**/
class Collect_common_factorsc:public map_function
{
    ex temex, temex2;
public:
    Collect_common_factorsc(){}
    ex operator()(const ex& _e);
    ~Collect_common_factorsc(){}
};

/**this Simplify2 function simplify only the algebraic expressions containing
fractional power of F,F_,Fd_,X_,Y_; such as: simplify (F^(1/3))^4 as F^(4/3).**/
ex Simplify2(const ex& expr_);

ex simplify(const ex& expr_, int rules = FuncSimp);
ex fullsimplify(const ex& expr_, int rules = FuncSimp);

class TrigArgSign_Complx:public map_function
{
    ex var_;
public:
    ex operator()(const ex & e);
    ~TrigArgSign_Complx(){}
};


/** expanding terms containing inverse power **/
class expandinv:public map_function
{
    exmap repls;
public:
    ex operator()(const ex& e);
    ~expandinv(){}
};

/**doing factors in power argument **/
class arguSimplify:public map_function
{
public:
    ex operator()(const ex& e);
    ~arguSimplify(){}
};



/** doing number simplify **/
class numSimplify:public map_function
{
public:
    map<numeric,numeric> primefactrs;
    ex getPrimefactors(const ex &e, const ex &fractimes);
    ex operator()(const ex& e);
    ~numSimplify(){}
};

/** basic simplification function to apply algebraic rules **/
//ex Simplify(const ex& _e, const int& rules = AlgSimp);

/** replacing the "pow" terms with created symbols, which have less degree than expandLevel and base is in add container. **/
class powBaseSubsLessThanDegLvl_1:public map_function
{
    unsigned j;
    ex expr;
    string str;
public:
    bool isNu;
    int addNum;
    exmap exprToSymMap;
    powBaseSubsLessThanDegLvl_1(){}
    void set(void)
    {
        j = 0;
        exprToSymMap.clear();addNum=0;isNu=false;
    }
    ex operator()(const ex& _e);
    ~powBaseSubsLessThanDegLvl_1(){}
};

/** replacing the "pow" terms with created symbols, which have less degree than expandLevel and base is in add container. **/
class fracPowBasSubsLvl_1:public map_function
{
    unsigned j;
    ex expr, tem;
    ex numer_denomClt;
    string str;
public:
    exmap baseCltLvl_1;
    fracPowBasSubsLvl_1(){}
    void set(void)
    {
        j = 0;
        baseCltLvl_1.clear();
    }
    ex operator()(const ex& e);
    ~fracPowBasSubsLvl_1(){}
};

/** replacing the "pow" terms with created symbols, which have less degree than expandLevel and base is in add container. **/
class fracPowBasSubsFactor:public map_function
{
    unsigned j;
    ex expr, tem;
    ex numer_denomClt;
    string str;
public:
    exmap baseClt;
    fracPowBasSubsFactor(){}
    void set(void)
    {
        j = 0;
        baseClt.clear();
    }
    ex operator()(const ex& e);
    ~fracPowBasSubsFactor(){}
};


/** replacing the "pow" terms with created symbols, which have less degree than expandLevel and base is in add container. **/
class powBaseSubsLessThanDeg:public map_function
{
    unsigned j;
    ex expr, tem;
    ex numer_denomClt;
    string str;
    powBaseSubsLessThanDegLvl_1 Lvl_1;
public:
    int addNum;
    exmap exprToSymMap;
    powBaseSubsLessThanDeg( unsigned j_): j(j_)
    {
        exprToSymMap.clear();
        addNum=0;
        Lvl_1.set();
    }
    ex operator()(const ex& _e);
    ~powBaseSubsLessThanDeg(){}
};


/** replacing base of fractional power with generated symbols  **/
class fracPowBasSubs:public map_function
{
    unsigned j;
    ex expr, tem;
    ex numer_denomClt;
    string str;
    fracPowBasSubsLvl_1 Lvl_1;
public:
    exmap baseClt;
    fracPowBasSubs(){}
    void set(void)
    {
        j = 0;
        Lvl_1.set();
        baseClt.clear();
    }
    ex operator()(const ex& e);
    ~fracPowBasSubs(){}
};

/** replacing some functions with generated symbols  **/
class funcSubs:public map_function
{
    unsigned j;
    ex expr,expr2,expr3;
    string str;
public:
    exmap baseClt;
    funcSubs(){j=0;baseClt.clear();}
    void set(void)
    {
        j = 0;
        baseClt.clear();
    }
    ex operator()(const ex& e);
    ~funcSubs(){}
};

/** Applying the simplification rules x^(3/2)=x*x^(1/2)  **/
class someMoreSimpRules:public map_function
{
    int iNum;

public:
    someMoreSimpRules(){}
    ex operator()(const ex& e);
    ~someMoreSimpRules(){}
};


extern simplifyc Simplify;
extern Collect_common_factorsc Collect_common_factors; // This collect all common factors including numerical numbers
extern numSimplify numSimplifye;
extern arguSimplify arguSimplifye;
extern expandinv expandinve;
extern fracPowBasSubs fracPowBasSubsE;
extern funcSubs funcSubsE;
extern someMoreSimpRules someMoreSimpRulesE;

#endif // SIMPLIFY_H_INCLUDED
