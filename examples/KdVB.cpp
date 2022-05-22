
/** @file NLS.cpp
 *
 *   This program solves KdVB equation
        {u_t} + u{u_x} - p{u_{xx}} + q{u_{xxx}} = 0,
    */



#include "GiNaCDE.h"


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"),p=reader("p"), q=reader("q"),
    k_0=reader("k_0"), k_1=reader("k_1"), A_0=reader("A_0"), A_1=reader("A_1"), 
    A_2=reader("A_2");  

    depend(u, {t,x});

    ex pde =Diff(u,t,1)+Diff(u,x,1)*u-p*Diff(u,x,2)+q*Diff(u,x,3);
    output = maple;  
    twcPhase = {lst{k_0,k_1},lst{}};
    degAcoeff = {2,1,1,A_2};
    ASolve=true;
    positivePart = true; 
    negativePart = true;
    paraInDiffSolve={};
    filename = "KdVB_mF.txt";
    desolve(pde, {u}, mF_expansion);

    return 0;
}
