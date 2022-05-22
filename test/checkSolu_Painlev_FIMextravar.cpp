
/** @file checkSolu_Painlev_FIMextravar.cpp
 *
 *   This program checks the solutions reported in Painlev_FIMextravar.txt file
    */



#include "GiNaCDE.h"


int main()
{

    // These declare dependent and independent variables of input diff. equ.
    const ex u=reader("u"),x=reader("x");

    // Make dependency of dependent variable on independent variables
    depend(u, {x});

    // Input equation is:
    const string DE = "Diff(u,x,2)+u^3*b+u*Diff(u,x,1)*a";


                    /************Check the solution #1***********/

    //Solutions of the algebraic equations are
    string algebraic_solutions = "{g_0==0,a_01==0,a_00==0,b==-1/2*(a+g_1)*g_1,a_02==1/2*a+1/2*g_1}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #1
    string solutions = "u = 2*(2*C_+(a+g_1)*x)^(-1)";

    //Following function, verify solution #1. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ex ret = checkSolu(DE,solutions,algebraic_solutions);
    if(ret!=_ex0)
    {
        return -1;
    }


                    /************Check the solution #4***********/

    //Solutions of the algebraic equations are
    algebraic_solutions = "{g_1==0,g_0==0,a_01==0,a_02==1/2*a,b==0}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #4
    solutions = "u = -tanh(C_+1/2*sqrt(-a_00*a)*sqrt(2)*x)^(-1)*(-a_00*a)^(-1/2)*sqrt(2)*a_00";

    //Following function, verify solution #4. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ret = checkSolu(DE,solutions,algebraic_solutions);
    if(ret!=_ex0)
    {
        return -1;
    }


    return 0;

}
