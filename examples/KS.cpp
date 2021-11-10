
/** @file KS.cpp
 *
 *   This program solves the Kudryashov<96>Sinelshchikov equation:
        {u_{3x}} + gu{u_x} - n{u_{2x}} - \left( {u{u_{2x}} + u_x^2} \right)d - k{u_x}{u_{2x}} - e\left( {u{u_{3x}} + {u_x}{u_{2x}}} \right) + {u_t} = 0,
    */



#include <GiNaCDE.h>


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"),
             A_0=reader("A_0"), A_1=reader("A_1"), A_2=reader("A_2"), A_3=reader("A_3"),
             g=reader("g"), n=reader("n"), d=reader("d"), e=reader("e"), k=reader("k");

    const ex pde = Diff(u,t,1)+g*u*Diff(u,x,1)+Diff(u,x,3)-e*Diff(u*Diff(u,x,2),x,1)-k*Diff(u,x,1)*Diff(u,x,2)-n*Diff(u,x,2)-d*Diff(u*Diff(u,x,1),x,1);

    depend(u, {t, x});

    //F-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{0,0}};
    degAcoeff=
    {3,A_0,A_1,A_2,A_3};
    ASolve=true;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="KS_Fexp.txt";
    desolve(pde, {u}, F_expansion);

    //mF-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{0,0}};
    degAcoeff=lst{2,A_0,A_1,A_2};
    ASolve=true;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="KS_mF.txt";
    desolve(pde, {u}, mF_expansion);

    //FIM method//
    twcPhase=lst{lst{k_0,k_1},lst{0,0}};
    paraInDiffSolve=lst{};
    filename="KS_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
