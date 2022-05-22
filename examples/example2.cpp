/** @file example2.cpp
 *
 *   Examples for GiNaCDE library. This program solve following pdes
        i.   Breaking Soliton equation.
        ii.  Kdv equation
        iii. Fisher equation
        iv.  Burger equation
        v.   Benjamin-Bona-Mahony equation
        vi.  Eckhaus equation
        vii. Boussines equation
        viii.Gardner equation
        ix.  Cahn-Allen equation
        x.   NLSE equation
        xi.  Kadomstev-Petviashvili equation
        xii. Duffing equation
        xiii.Landau-Ginzburg-Higgs equation
        xiv. Kudryashov-Sinelshchikov  equation
        xv.  Generalized Camassa-Holm equation
        xvi. Nonlinear Telegraph Equation      */


#include "GiNaCDE.h"


int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), y=reader("y"),z=reader("z"),a=reader("a"), b=reader("b"),
             c=reader("c"),p=reader("p"),q=reader("q"),k=reader("k"),w=reader("w"),k_0=reader("k_0"),k_1=reader("k_1"),
             k_2=reader("k_2"),k_3=reader("k_3"),A_0=reader("A_0"),A_1=reader("A_1"),A_2=reader("A_2"),A_3=reader("A_3"),
             A_4=reader("A_4"),A_5=reader("A_5"),A_6=reader("A_6");
    ex pde;

    depend(u, {t, x, y, z});

    pde = Diff(Diff(u,x,1),t,1) - 4*Diff(u,x,1)*Diff(Diff(u,x,1),y,1) - 2*Diff(u,y,1)*Diff(u,x,2) + Diff(Diff(u,x,3),y,1);  // Breaking Soliton equation
    output = maple;  // Results are saved in maple language
    twcPhase=lst{lst{k_0,k_1,k_2},lst{}};
    degAcoeff = lst{2,0,A_1,A_2};
    filename = "breakingSoliton_mF.txt";
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,1)-u*u*a*Diff(u,x,1)+Diff(u,x,3);  // KDV equation
    output = mathematica; // Results are saved in mathematica language
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "KDV_FIM1.txt";
    desolve(pde,{u}, FIM);
    NValue = 2;  // defined value of N
    output = maple;
    filename = "KDV_FIM2.txt";
    desolve(pde, {u}, FIM);
    output = mathematica;
    filename = "KDV_Fex.txt";
    paraInDiffSolve = lst{a};
    degAcoeff = lst{4,A_0,0,A_2,0,A_4};
    desolve(pde, {u}, F_expansion);

    pde = Diff(u,t,1)-Diff(u,x,2)-u+pow(u,2); // Fisher equation
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "fisher_Fex.txt";
    ASolve = true;
    degAcoeff = lst{4,A_0,0,A_2,0,A_4};
    desolve(pde, {u}, F_expansion);

    pde = Diff(u,t,1)+Diff(u,x,1)*u-a*Diff(u,x,2); // Burger equation
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "burger_mF.txt";
    degAcoeff = lst{2,A_0,A_1,A_2};
    ASolve = true;
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,1)+Diff(u,x,1)+u*Diff(u,x,1)-Diff(Diff(u,x,2),t,1);  // Benjamin-Bona-Mahony equation
    NValue = 2;
    output = mathematica;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "benjamin_FIM.txt";
    desolve(pde, {u}, FIM);
    filename = "benjamin_Fex.txt";
    degAcoeff = lst{4,0,0,A_2,0,A_4};
    ASolve = true;
    positivePart = true;
    negativePart = false;
    desolve(pde, {u}, F_expansion);

    pde = I*Diff(u,t,1) + Diff(u,x,2) + 2*u*Diff(u*conjugate(u),x,1) + u*u*conjugate(u)*conjugate(u)*u;  // Eckhaus equation
    output = mathematica;
    twcPhase=lst{lst{-2*k*a,k},lst{b,a}};
    filename = "eckhaus_FIM.txt";
    desolve(pde, {u}, FIM);
    output = maple;
    degAcoeff = lst{2,A_0,0,A_2};
    filename = "eckhaus_mF.txt";
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,2)-Diff(u,x,2)-3*Diff(u*u,x,2)-Diff(u,x,4);  // Boussines equation
    output = maple;
    twcPhase = lst{lst{k_0,k_1},lst{}};
    degAcoeff = lst{2,A_0,A_1,A_2};
    ASolve = true;
    filename = "boussines_mF.txt";
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,1)+2*a*u*Diff(u,x,1)-3*b*u*u*Diff(u,x,1)+Diff(u,x,3);  // Gardner equation
    output = mathematica;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "gardner_Fex.txt";
    ASolve = true;
    degAcoeff = lst{4,0,0,A_2,A_3,A_4};
    desolve(pde, {u}, F_expansion);
    filename = "gardner_FIM.txt";
    NValue = 2;
    desolve(pde, {u}, FIM);

    pde = Diff(u,t,1) - Diff(u,x,2) + u*u*u - u;  // Cahn-Allen equation
    output = ginac;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "cahnAllen_mF.txt";
    degAcoeff = lst{2,0,A_1,A_2};
    positivePart = true;
    negativePart = false;
    desolve(pde, {u}, mF_expansion);
    output = ginac;
    filename = "cahnAllen_FIM.txt";
    desolve(pde, {u}, FIM);

    pde = I*Diff(u,t,1) + p*Diff(u,x,2) + q*u*u*conjugate(u);  // NLSE equation
    output = maple;
    twcPhase=lst{lst{-2*p*a,1},lst{b,a}};
    filename = "NLSE_Fex.txt";
    degAcoeff = lst{4,A_0,0,A_2,0,A_4};
    desolve(pde, {u}, F_expansion);
    filename = "NLSE_FIM.txt";
    desolve(pde, {u}, FIM);

    pde = Diff(Diff(u,t,1),x,1)+6*Diff(u,x,1)*Diff(u,x,1)+6*u*Diff(u,x,2)+Diff(u,x,3)-Diff(u,y,2)-Diff(u,z,2);  // Kadomstev-Petviashvili equation
    output = mathematica;
    twcPhase=lst{lst{k_0,k_1,k_2,k_3},lst{}};
    degAcoeff = lst{2,0,A_1,A_2};
    filename = "Kadomstev_mF.txt";
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,2)+b*u+c*pow(u,3);  // Duffing equation
    output = mathematica;
    twcPhase=lst{lst{},lst{}};
    filename = "Duffing_FIM.txt";
    paraInDiffSolve = lst{b,c};
    desolve(pde, {u}, FIM);
    NValue = 2;
    filename = "Duffing_FIM2.txt";
    paraInDiffSolve = lst{b,c};
    desolve(pde, {u}, FIM);

    pde = Diff(u,t,2)-Diff(u,x,2)+pow(b,2)*u+pow(c,2)*pow(u,3);  // Landau-Ginzburg-Higgs equation
    output = maple;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "Landau-Ginzburg-Higgs_FIM.txt";
    paraInDiffSolve = lst{b,c};
    desolve(pde, {u}, FIM);
    NValue=2;
    filename = "Landau-Ginzburg-Higgs_FIM2.txt";
    paraInDiffSolve = lst{b,c};
    desolve(pde, {u}, FIM);

    pde = Diff(u,t,1)+u*Diff(u,x,1)+Diff(u,x,3)-Diff(u*Diff(u,x,2),x,1)-Diff(u,x,1)*Diff(u,x,2)-Diff(u,x,2)-Diff(u*Diff(u,x,1),x,1);  // Kudryashov-Sinelshchikov  equation
    output = mathematica;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff = lst{2,0,A_1,A_2};
    NValue = 2;
    filename = "Kudryashov-Sinelshchikov_mF.txt";
    desolve(pde, {u}, mF_expansion);

    pde = Diff(u,t,1)+2*k*Diff(u,x,1)-Diff(Diff(u,x,2),t,1)+a*u*Diff(u,x,1)-2*Diff(u,x,1)*Diff(u,x,2)-u*Diff(u,x,3);  // Generalized Camassa-Holm equation
    output = maple;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff = lst{3,0,A_1,0,A_3};
    NValue = 2;
    filename = "Generalized_Camassa-Holm_mF.txt";
    paraInDiffSolve = lst{k,a};
    desolve(pde, {u}, mF_expansion);


    pde = Diff(u,t,2)-Diff(u,x,2)+Diff(u,t,1)+a*u+b*pow(u,3);  //  Nonlinear Telegraph Equation
    output=maple;
    twcPhase=lst{lst{k_0,k_1},lst{0,0}};
    degAcoeff=
    {3,A_0,A_1,A_2,A_3};
    ASolve=false;
    positivePart=true;
    negativePart=true;
    paraInDiffSolve=lst{};
    filename="telegraph_Fexp.txt";
    desolve(pde, {u}, F_expansion);
    NValue=1;
    filename="telegraph_FIM.txt";
    desolve(pde, {u}, FIM);

    return 0;

}
