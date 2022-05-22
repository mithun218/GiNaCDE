
/** @file kerrNLS.cpp
 *
 *   This program solves the perturbed NLS equation with Kerr law nonlinearity:
        I{u_t} + {u_{2x}} + A|u{|^2}u + I\left( {{G_1}{u_{3x}} + {G_2}|u{|^2}{u_x} + {G_3}{{\left( {|u{|^2}} \right)}_x}u} \right) = 0,
    */



#include "GiNaCDE.h"


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), k_0=reader("k_0"), k_1=reader("k_1"), p_0=reader("p_0"), p_1=reader("p_1"),
             A_0=reader("A_0"), A_1=reader("A_1"), A_2=reader("A_2"), A_3=reader("A_3"), A_4=reader("A_4"),
             A=reader("A"), G_1=reader("G_1"), G_2=reader("G_2"), G_3=reader("G_3");

    const ex pde = I*Diff(u,t,1)+Diff(u,x,2)+A*conjugate(u)*pow(u,2)+I*(G_1*Diff(u,x,3)+G_2*conjugate(u)*u*Diff(u,x,1)+G_3*Diff(u*conjugate(u),x,1)*u);

    depend(u, {t, x});

    //F-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{4,0,0,A_2,A_3,A_4};
    ASolve=true;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="kerrNLS_Fexp.txt";
    desolve(pde, {u}, F_expansion);

    //mF-expansion method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    degAcoeff=lst{2,A_0,A_1,A_2};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="kerrNLS_mF.txt";
    desolve(pde, {u}, mF_expansion);

    //FIM method//
    twcPhase=lst{lst{k_0,k_1},lst{p_0,p_1}};
    paraInDiffSolve=lst{};
    filename="kerrNLS_FIM.txt";
    desolve(pde, {u}, FIM);


    return 0;

}
