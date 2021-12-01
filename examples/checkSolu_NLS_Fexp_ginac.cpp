
/** @file checkSolu_NLS_Fexp_ginac.cpp
 *
 *   This program checks the solutions reported in NLS_Fexp_ginac.txt file
    */



#include <GiNaCDE/GiNaCDE.h>


int main()
{

    // These declare dependent and independent variables of input diff. equ.
    const ex u=reader("u"), t=reader("t"), x=reader("x");

    // Make dependency of dependent variable on independent variables
    depend(u, {t, x});

    // Input equation is:
    const string DE = "I*Diff(u,t,1)-Diff(u,x,2)*p+u^2*q*conjugate(u)";

    //The input NLPDE is complex.
    //We derive solutions of real part of Diff. Equ. with the condition:
    string conditions = "p_1 = 1/2*p^(-1)*k_1^(-1)*k_0";


                    /************Check the solution #2***********/

    //Solutions of the algebraic equations are
    string algebraic_solutions = "{a_0==0,p_0==-1/4*(4*p^2*k_1^4*A_2-k_0^2)*p^(-1)*k_1^(-2),b_1==-q^(-1/2)*sqrt(2)*sqrt(p)*k_1*sqrt(A_0),a_1==0}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #2
    string solutions = "u=2*(exp(2*C_*sqrt(A_2))*A_0-exp(2*(k_1*x+k_0*t)*sqrt(A_2)))^(-1)*q^(-1/2)*sqrt(2)*sqrt(p)*k_1*sqrt(A_2)*exp(-(C_+k_1*x+k_0*t)*sqrt(A_2))^(-1)*sqrt(A_0)*exp(-(1/4*I)*(4*p^2*k_1^4*A_2-k_0^2)*p^(-1)*k_1^(-2)*t+(1/2*I)*p^(-1)*k_1^(-1)*k_0*x)";

    //Following function, verify solution #2. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ex ret = checkSolu(DE,solutions,algebraic_solutions,conditions);
    if(ret==_ex0)
    {
        cout<<"solution #2 is correct:)"<<endl;
        cout<<ret<<endl;
    }
    else
    {
        cout<<"solution #2 is not correct."<<endl;
        cout<<ret<<endl;
    }

                    /************Check the solution #13***********/

    //Solutions of the algebraic equations are
    algebraic_solutions = "{a_0==0,b_1==q^(-1/2)*sqrt(2)*sqrt(p)*k_1*sqrt(A_0),p_0==-1/4*(4*p^2*k_1^4*A_2-k_0^2)*p^(-1)*k_1^(-2),a_1==0}";

    //solution(s) of input Diff. Equ. is (are)=>
    //solution #13
    solutions = "u = -sinh(C_+(k_1*x+k_0*t)*sqrt(A_2))^(-1)*q^(-1/2)*(A_2*A_0)^(-1/2)*sqrt(2)*sqrt(p)*k_1*A_2*sqrt(A_0)*exp(-(1/4*I)*(4*p^2*k_1^4*A_2-k_0^2)*p^(-1)*k_1^(-2)*t+(1/2*I)*p^(-1)*k_1^(-1)*k_0*x)";

    //Following function, verify solution #13. If the following function returns 0, the solution is correct.
    //Otherwise, the solution is not valid.
    ret = checkSolu(DE,solutions,algebraic_solutions,conditions);
    if(ret==_ex0)
    {
    cout<<"solution #13 is correct:)"<<endl;
    cout<<ret<<endl;
    }
    else
    {
    cout<<"solution #13 is not correct."<<endl;
    cout<<ret<<endl;
    }


    return 0;

}
