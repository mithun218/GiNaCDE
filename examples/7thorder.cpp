
/** @file 7thorder.cpp
 *
 *   This program solves Seventh-order Sawada-Kotera equation (sSK):
        {u_t} + {(63{u^4} + 63(2{u^2}{u_{2x}} + uu_x^2) + 21(u{u_{4x}} + u_{2x}^2 + {u_x}{u_{3x}}) + {u_{6x}})_x} = 0
    */



#include <GiNaCDE/GiNaCDE.h>


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             A_0=reader("A_0"), A_2=reader("A_2"), A_4=reader("A_4");

    const ex pde = Diff(u,t,1)+Diff((63*pow(u,4)+63*(2*pow(u,2)*Diff(u,x,2)+u*pow(Diff(u,x,1),2)))+21*(u*Diff(u,x,4)+pow(Diff(u,x,2),2)+Diff(u,x,1)*Diff(u,x,3))+Diff(u,x,6),x,1);

    depend(u, {t, x});

    //F-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff=lst{4,A_0,0,A_2,0,A_4};
    ASolve=true;
    positivePart=true;
    negativePart=false;
    paraInDiffSolve=lst{};
    filename="7thorder_Fexp.txt";
    desolve(pde, {u}, F_expansion);

    //mF-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff=lst{2,0,1,1};
    ASolve=false;
    positivePart=true;
    negativePart=false;
    paraInDiffSolve=lst{};
    filename="7thorder_mF.txt";
    desolve(pde, {u}, mF_expansion);


    //FIM method//
    twcPhase=lst{lst{k_0,k_1},lst{}};
    paraInDiffSolve=lst{};
    filename="7thorder_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
