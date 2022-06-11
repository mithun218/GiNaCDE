
/** @file test2.cpp
 *
 *   Program to test GiNaCDE library. This program solve following pdes
        i.  Cahn-Allen equation
        ii. Duffing equation
        iii.Nonlinear Telegraph Equation 
        iv. Kudryashov-Sinelshchikov  equation 
        v.  Benjamin-Bona-Mahony equation    */



#include "GiNaCDE.h"
//#include <conio.h>

int main()
{
    const ex u=reader("u"), t=reader("t"), x=reader("x"), y=reader("y"),z=reader("z"),a=reader("a"), b=reader("b"),
             c=reader("c"),p=reader("p"),q=reader("q"),k=reader("k"),w=reader("w"),k_0=reader("k_0"),k_1=reader("k_1"),
             k_2=reader("k_2"),k_3=reader("k_3"),A_0=reader("A_0"),A_1=reader("A_1"),A_2=reader("A_2"),A_3=reader("A_3"),
             A_4=reader("A_4"),A_5=reader("A_5"),A_6=reader("A_6");
    ex pde, res;
    size_t solu_num1,solu_num2;
    stringstream diffStr,algSoluStr,diffSoluStr;
    string str;

    depend(u, {t, x, y, z});


    pde = Diff(u,t,1) - Diff(u,x,2) + u*u*u - u;  // Cahn-Allen equation
    output = ginac;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    filename = "cahnAllen_mF.txt";
    degAcoeff = lst{2,0,A_1,A_2};
    positivePart = true;
    negativePart = false;
    desolve(pde, {u}, mF_expansion,true);
    /* Checking all solutions*/
     diffStr.str("");
     diffStr<<pde;
     solu_num1 = (solutionClt).size();
     for(size_t i=0;i<solu_num1;i++)
     {
         algSoluStr.str("");
         algSoluStr<<solutionClt[i][0];
         solu_num2 = nops(solutionClt[i]);
         for(size_t j=1;j<solu_num2;j++)
         {
             diffSoluStr.str("");
             diffSoluStr<<solutionClt[i][j];
             res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
             if(res!=_ex0)
                 return -1;
         }

     }


    pde = Diff(u,t,2)+b*u+c*pow(u,3);  // Duffing equation
    output = ginac;
    twcPhase=lst{lst{},lst{}};
    filename = "Duffing_FIM.txt";
    paraInDiffSolve = lst{b,c};
    desolve(pde, {u}, FIM,true);
    /* Checking all solutions*/
     diffStr.str("");
     diffStr<<pde;
     solu_num1 = (solutionClt).size();
     for(size_t i=0;i<solu_num1;i++)
     {
         algSoluStr.str("");
         algSoluStr<<solutionClt[i][0];
         solu_num2 = nops(solutionClt[i]);
         for(size_t j=1;j<solu_num2;j++)
         {
             diffSoluStr.str("");
             diffSoluStr<<solutionClt[i][j];
             res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
             if(res!=_ex0)
                 return -1;
         }

     }


    pde = Diff(u,t,2)-Diff(u,x,2)+Diff(u,t,1)+a*u+b*pow(u,3);  //  Nonlinear Telegraph Equation
    output=ginac;
    twcPhase=lst{lst{k_0,k_1},lst{0,0}};
    NValue=1;
    filename="telegraph_FIM.txt";
    desolve(pde, {u}, FIM,true);
    /* Checking all solutions*/
     diffStr.str("");
     diffStr<<pde;
     solu_num1 = (solutionClt).size();
     for(size_t i=0;i<solu_num1;i++)
     {
         algSoluStr.str("");
         algSoluStr<<solutionClt[i][0];
         solu_num2 = nops(solutionClt[i]);
         for(size_t j=1;j<solu_num2;j++)
         {
             diffSoluStr.str("");
             diffSoluStr<<solutionClt[i][j];
             res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
             if(res!=_ex0)
                 return -1;
         }

     }


   pde = Diff(u,t,1)+u*Diff(u,x,1)+Diff(u,x,3)-Diff(u*Diff(u,x,2),x,1)-Diff(u,x,1)*Diff(u,x,2)-Diff(u,x,2)-Diff(u*Diff(u,x,1),x,1);  // Kudryashov-Sinelshchikov  equation
    output = ginac;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff = lst{2,0,A_1,A_2};
    NValue = 2;
    filename = "Kudryashov-Sinelshchikov_mF.txt";
    desolve(pde, {u}, mF_expansion,true);
    /* Checking all solutions*/
     diffStr.str("");
     diffStr<<pde;
     solu_num1 = (solutionClt).size();
     for(size_t i=0;i<solu_num1;i++)
     {
         algSoluStr.str("");
         algSoluStr<<solutionClt[i][0];
         solu_num2 = nops(solutionClt[i]);
         for(size_t j=1;j<solu_num2;j++)
         {
             diffSoluStr.str("");
             diffSoluStr<<solutionClt[i][j];
             res = checkSolu(diffStr.str(),diffSoluStr.str(),algSoluStr.str());
             if(res!=_ex0)
                 return -1;
         }

     }


    pde = Diff(u,t,1)+Diff(u,x,1)+u*Diff(u,x,1)-Diff(Diff(u,x,2),t,1);  // Benjamin-Bona-Mahony equation
    output = ginac;
    twcPhase=lst{lst{k_0,k_1},lst{}};
    degAcoeff = lst{4,0,0,A_2,0,A_4};
    ASolve = true;
    positivePart = true;
    negativePart = false;
    filename = "benjamin_Fex.txt";
    desolve(pde, {u}, F_expansion,true);
    /* Checking all solutions*/
     diffStr.str("");
     diffStr<<pde;
     solu_num1 = (solutionClt).size();
     for(size_t i=0;i<solu_num1;i++)
     {
         algSoluStr.str("");
         algSoluStr<<solutionClt[i][0];
         solu_num2 = nops(solutionClt[i]);
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
