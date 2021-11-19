
/** @file checkSolu_cahnAllen_mF.cpp
 *
 *   This program checks the solutions reported in cahnAllen_mF.txt file
    */



#include <GiNaCDE/GiNaCDE.h>


int main()
{

    // These declare dependent and independent variables of input diff. equ.
    const ex u=reader("u"), t=reader("t"), x=reader("x");

    // Make dependency of dependent variable on independent variables
    depend(u, {t, x});

    // Input equation is:
    const string DE = "-u-Diff(u,x,2)+Diff(u,t,1)+u^3";


                        /************Check the solution #1***********/

    //Solutions of the algebraic equations are
    string algebraic_solutions = "{k_0==-3/2*A_1^(-1),a_0==-1,a_1==-A_2*A_1^(-1),k_1==1/2*sqrt(2)*A_1^(-1)}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #1
    string solutions = "u = -1+(cosh(A_1*C_+1/2*(sqrt(2)*A_1^(-1)*x-3*A_1^(-1)*t)*A_1)+sinh(A_1*C_+1/2*(sqrt(2)*A_1^(-1)*x-3*A_1^(-1)*t)*A_1))*(-1+cosh(A_1*C_+1/2*(sqrt(2)*A_1^(-1)*x-3*A_1^(-1)*t)*A_1)*A_2+sinh(A_1*C_+1/2*(sqrt(2)*A_1^(-1)*x-3*A_1^(-1)*t)*A_1)*A_2)^(-1)*A_2";

    //Following function, verify solution #1. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ex ret = checkSolu(DE,solutions,algebraic_solutions);
    if(ret==_ex0)
    {
        cout<<"solution #1 is correct:)"<<endl;
        cout<<ret<<endl;
    }
    else
    {
        cout<<"solution #1 is not correct."<<ret<<endl;
        cout<<ret<<endl;
    }

                        /************Check the solution #20***********/

    //Solutions of the algebraic equations are
    algebraic_solutions = "{k_0==3/2*A_1^(-1),a_1==A_2*A_1^(-1),a_0==0,k_1==-1/2*sqrt(2)*A_1^(-1)}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #20
    solutions = "u = -1/2*(tanh(-1/4*(sqrt(2)*A_1^(-1)*x-3*A_1^(-1)*t)*A_1+C_)*A_2^(-1)*A_1+A_2^(-1)*A_1)*A_2*A_1^(-1)";

    //Following function, verify solution #20. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ret = checkSolu(DE,solutions,algebraic_solutions);
    if(ret==_ex0)
    {
        cout<<"solution #20 is correct:)"<<endl;
        cout<<ret<<endl;
    }
    else
    {
        cout<<"solution #20 is not correct."<<endl;
        cout<<ret<<endl;
    }


    return 0;

}
