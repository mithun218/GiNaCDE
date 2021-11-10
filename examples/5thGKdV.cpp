
/** @file 5thGKdV.cpp
 *
 *   This program solves fifth-order KdV equation (gKdV):
        pu{u_{3x}} + q{u_x}{u_{2x}} + r{u^2}{u_x} + {u_{5x}} + {u_t} = 0,
    */



#include <GiNaCDE.h>


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             A_1=reader("A_1"), A_2=reader("A_2"), A_3=reader("A_3"), p=reader("p"), q=reader("q"), r=reader("r");

    const ex pde = Diff(u,t,1)+p*u*Diff(u,x,3)+q*Diff(u,x,1)*Diff(u,x,2)+r*pow(u,2)*Diff(u,x,1)+Diff(u,x,5);

    depend(u, {t, x});

    //F-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff=lst{3,0,A_1,A_2,A_3};
    ASolve=false;
    positivePart=true;
    negativePart=false;
    paraInDiffSolve=lst{r};
    filename="5thGKdV_Fexp.txt";
    desolve(pde, {u}, F_expansion);

    //mF-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff=lst{2,1,1,1};
    ASolve=false;
    positivePart=true;
    negativePart=false;
    paraInDiffSolve=lst{q};
    filename="5thGKdV_mF.txt";
    desolve(pde, {u}, mF_expansion);

    //FIM method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    paraInDiffSolve=lst{};
    filename="5thGKdV_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
