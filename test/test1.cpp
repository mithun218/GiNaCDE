
/** @file test1.cpp
 *
 *   Program to test GiNaCDE library. This solve following odes
        i.   linear ode - the damped harmonic oscillator.
        ii.  nonlinear ode - modified Painlev-Ince equation */




#include "GiNaCDE.h"

 int main()
 {
    const ex u=reader("u"),x=reader("x"),a=reader("a"),b=reader("b"),
             c=reader("c"),w=reader("w"),A_0=reader("A_0"),A_1=reader("A_1"),A_2=reader("A_2");
    ex ode,res;
    size_t solu_num1,solu_num2;
    stringstream diffStr,algSoluStr,diffSoluStr;
    string str;

    depend(u, {x});

    twcPhase = lst{lst{},lst{}};

    ode = Diff(u,x,2) + Diff(u,x,1) + w*w*u; // the damped harmonic oscillator
    output = maple;// Outputs are saved in maple format;
    filename = "damped_FIM.txt";
    desolve(ode,{u},FIM,true);
    /* Checking all solutions*/
    diffStr.str("");
    diffStr<<ode;
    solu_num1 = sizeof(solutionClt);
    for(size_t i=1;i<solu_num1;i++)
    {
        algSoluStr.str("");
        algSoluStr<<solutionClt[i][0];
        solu_num2 = sizeof(solutionClt[i]);
        for(size_t j=1;j<solu_num2;j++)
        {
            diffSoluStr.str("");
            diffSoluStr<<solutionClt[i][j];
            res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
            if(res!=_ex0)
                return -1;
        }

    }



    ode = Diff(u,x,2) + a*u*Diff(u,x,1) + b*u*u*u; // modified Painlev-Ince equation
    filename = "Painlev_FIM.txt";
    desolve(ode,{u},FIM,true);
    output = ginac;
    filename = "Painlev_FIMextravar.txt";
    paraInDiffSolve = lst{a,b};
    desolve(ode,{u},FIM,true);
    /* Checking all solutions*/
    diffStr.str("");
    diffStr<<ode;
    solu_num1 = sizeof(solutionClt);
    for(size_t i=1;i<solu_num1;i++)
    {
        algSoluStr.str("");
        algSoluStr<<solutionClt[i][0];
        solu_num2 = sizeof(solutionClt[i]);
        for(size_t j=1;j<solu_num2;j++)
        {
            diffSoluStr.str("");
            diffSoluStr<<solutionClt[i][j];
            res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
            if(res!=_ex0)
                return -1;
        }

    }

    return 0;

 }
