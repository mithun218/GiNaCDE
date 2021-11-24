
/** @file Eckhaus.cpp
 *
 *   This program solves Eckhaus equation:
        I{u_t} + {u_{xx}} + 2{\left( {{{\left| u \right|}^2}} \right)_x}u + {\left| u \right|^4}u = 0
    */



#include <GiNaCDE/GiNaCDE.h>


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             p_0=reader("p_0"), p_1=reader("p_1"), A_0=reader("A_0"), A_1=reader("A_1"),
             A_2=reader("A_2");

    const ex pde = I*Diff(u,t,1) + Diff(u,x,2) + 2*u*Diff(u*conjugate(u),x,1) + u*u*conjugate(u)*conjugate(u)*u;

    depend(u, {t, x});

    //F-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{1,A_0,A_1};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="EckhausFexp.txt";
    desolve(pde, {u}, F_expansion);

    //mF-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{2,0,A_1,A_2};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="EckhausmF.txt";
    desolve(pde, {u}, mF_expansion);

    //FIM method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    paraInDiffSolve=lst{};
    filename="EckhausFIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
