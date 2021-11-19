
/** @file NLS.cpp
 *
 *   This program solves one dimensional cubic nonlinear Schr\"odinger (NLS) equation:
        Iu_t-pu_{xx}+q{|u|}^2u=0,
    */



#include <GiNaCDE.h>


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             p_0=reader("p_0"), p_1=reader("p_1"), A_0=reader("A_0"),A_2=reader("A_2"),
             p=reader("p"),q=reader("q");

    const ex pde = I*Diff(u,t,1) - p*Diff(u,x,2) + q*u*u*conjugate(u);

    depend(u, {t, x});

    output=maple;
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{2,A_0,0,A_2};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="NLS_Fexp(maple).txt";
    desolve(pde,{u},F_expansion);
    output=ginac;
    filename="NLS_Fexp(ginac).txt";
    desolve(pde,{u},F_expansion);

    output=mathematica;
    filename="NLS_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
