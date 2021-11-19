
/** @file test1.cpp
 *
 *   Program to test GiNaCDE library. This solve following odes
        i.   linear ode - the damped harmonic oscillator.
        ii.  nonlinear ode - modified Painlev-Ince equation
        iii. Kudryashov-Sinelshchikov equation */




#include "GiNaCDE.h"

 int main()
 {
    const ex u=reader("u"),x=reader("x"),a=reader("a"),b=reader("b"),
             c=reader("c"),w=reader("w"),A_0=reader("A_0"),A_1=reader("A_1"),A_2=reader("A_2");
    ex ode;

    depend(u, {x});

    twcPhase = lst{lst{},lst{}};

    ode = Diff(u,x,2) + Diff(u,x,1) + w*w*u; // the damped harmonic oscillator
    output = maple;// Outputs are saved in maple format;
    filename = "damped_FIM.txt";
    desolve(ode,{u},FIM,true);


    ode = Diff(u,x,2) + a*u*Diff(u,x,1) + b*u*u*u; // modified Painlev-Ince equation
    filename = "Painlev_FIM.txt";
    desolve(ode,{u},FIM,true);
    output = ginac;
    filename = "Painlev_FIMextravar.txt";
    paraInDiffSolve = lst{a,b};
    desolve(ode,{u},FIM,true);


    ode = c*Diff(u,x,1)-a*(1-u)*Diff(u,x,1)-Diff(u,x,3)+Diff((1-u)*Diff(u,x,2),x,1)-b*Diff(u,x,1)*Diff(u,x,2);  // Kudryashov-Sinelshchikov equation
    output = maple;
    NValue = 1;
    degAcoeff = lst{2,1,1,A_2};
    filename = "kudryashov_mF.txt";
    paraInDiffSolve = lst{a,b,c};
    desolve(ode, {u}, mF_expansion,true);



    return 0;

 }
