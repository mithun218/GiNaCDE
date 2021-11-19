
/** @file gtools.cpp
 *
 *  A tool for solving differential equations using GiNaCDE library.*/

#include "GiNaCDE.h"

using namespace std;
using namespace GiNaC;

int main()
{
    cout << "gtools - a tool for solving differential equation (GiNaCDE V" << GINACDE_VERSION_MAJOR << "." << GINACDE_VERSION_MINOR << "." <<GINACDE_VERSION_MAINTENANCE<<")"<< endl<< endl;

    string inputstr;
    int method;
    ex inputdiff, dpndt_var;
    lst extravar, indpndt_vars;
    parser reader;

    cout << "Type h for help." << endl << endl;


    do
    {
        cout << "Input dependent variable: ";
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Please provide dependent variable of input diff. equ." ;
        }
        cout << endl;
    }while(inputstr == "h");
    dpndt_var = reader(inputstr);


    do
    {
        cout << "Input independent variables: ";
        cin >> inputstr;
        if(inputstr == "h" || !is_a<lst>(reader(inputstr)))
        {
            cout << endl << "Please provide independent variables of input diff. equ. between curly\n"
                            "brackets. Example: {t,x,y,z}." ;
        }
        cout << endl;
    }while(inputstr == "h" || !is_a<lst>(reader(inputstr)));
    indpndt_vars = ex_to<lst>(reader(inputstr));
    for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
    {
        depend(dpndt_var, {*it});
    }

    do
    {
        cout << "Input differential equation: ";
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Example of an input diff: Diff(u,t,1)+Diff(u,x,1)+u*Diff(u,x,1)-Diff(Diff(u,x,2),t,1)" ;
        }
        cout << endl;
    }while(inputstr == "h");

    inputdiff = reader(inputstr);

    do
    {
        cout << "Output format for saving results: " ;
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Type m for output format in maple. \n"
                            "Type M for output format in mathematica. \n"
                            "Type g for output format in ginac.";
        }
        cout << endl;
    }while(inputstr != "m" && inputstr != "maple" && inputstr != "M" && inputstr != "mathematica");
    if(inputstr == "m" || inputstr == "maple")
        output = maple;
    else if(inputstr == "M" || inputstr == "mathematica")
        output = mathematica;
    else
        output = ginac;

    do
    {
        cout << "Do you assign value of N?: " ;
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Type y for yes.\n"
                            "Type n for no.\n"
                            "For FIM, 'N' represents number of terms in sum(a_i(X)*Y^i,i=0..N) where X=u(x), Y=diff(u(x),x). For FIM, N = 1 and N = 2 are only allowed. "
                            "For F-expansion and modified F-expansion methods, the solutions of input NLPDE is expressed by a finite power series where presents 'N+1' terms.\n"
                            "By default, for FIM, N=1 and for F-expansion, modified F-expansion methods 'N' is evaluated automatically if N is not assigned to any value.";
        }
        cout << endl;
    }while(inputstr != "y" && inputstr != "n");
    if(inputstr == "y")
    {
        do
        {
            cout << "N: ";
            cin >> inputstr;
        }while(!has_only_digits(inputstr));
        stringstream ss(inputstr);
        int temNValue;
        ss >> temNValue;
        NValue = temNValue;
    }

    twcPhase = {lst{},lst{}};

    do
    {
        cout << "Provide constants in the traveling wave coordinate: " ;
        cin >> inputstr;
        if(inputstr == "h" || !is_a<lst>(reader(inputstr)))
        {
            cout << endl << "ex: {k_0,k_1}";
        }
        cout << endl;
    }while(inputstr == "h" || !is_a<lst>(reader(inputstr)));
    twcPhase[0]=(ex_to<lst>(reader(inputstr)));

    conjuFree conjuFree;
    replaceI replaceI;
    const ex replaceIex = replaceI(inputdiff),
             conjuFreeex = conjuFree(inputdiff);
    if(replaceIex != inputdiff || conjuFreeex != inputdiff)
    {
        do
        {
            cout << "Provide constants in phase part: " ;
            cin >> inputstr;
            if(inputstr == "h" || !is_a<lst>(reader(inputstr)))
            {
                cout << endl << "ex: {kp_0,kp_1}";
            }
            cout << endl;
        }while(inputstr == "h" || !is_a<lst>(reader(inputstr)));
        twcPhase[1]=(ex_to<lst>(reader(inputstr)));
    }


    do
    {
        cout << "Methods for solving differential equation : " ;
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Type f for first integral method.\n"
                            "Type F for F-expansion method.\n"
                            "Type mF for modified F-expansion method." ;
        }
        cout << endl;
    }while(inputstr != "f" && inputstr != "FIM" && inputstr != "F" && inputstr != "F_expansion" && inputstr != "mF" && inputstr != "mF_expansion" );
    if(inputstr == "f" || inputstr == "FIM")
        method = FIM;
    else if(inputstr == "F" || inputstr == "F_expansion")
    {
        method = F_expansion;
    }
    else if(inputstr == "mF" || inputstr == "mF_expansion")
    {
        method = mF_expansion;
    }

    if(method == F_expansion || method == mF_expansion)
    {
        do
        {
            cout << "Provide highest integer delta and coefficients in first-order nonlinear ode (A.E.): " ;
            cin >> inputstr;
            if(inputstr == "h" || !is_a<lst>(reader(inputstr)))
            {
                cout << endl << "ex: {2,A_0,A_1,A_2}";
            }
            cout << endl;
        }while(inputstr == "h" || !is_a<lst>(reader(inputstr)));
        degAcoeff = ex_to<lst>(reader(inputstr));

        do
        {
            cout << "ASolve: ";
            cin >> inputstr;
            if(inputstr == "h")
            {
                cout << endl <<"Type y for solving algebraic system for parameters A_i (i=0,1,..delta) in A.E.., "
                               "otherwise type n." ;
            }
            cout << endl;
        }while(inputstr != "y" && inputstr != "n");

        if(inputstr == "y")
            ASolve = true;
        else
            ASolve = false;


        do
        {
            cout << "positivePart: ";
            cin >> inputstr;
            if(inputstr == "h")
            {
                cout << endl <<"Type y for retaining positive part in solutions,, "
                               "otherwise type n." ;
            }
            cout << endl;
        }while(inputstr != "y" && inputstr != "n");

        if(inputstr == "y")
            positivePart = true;
        else
            positivePart = false;


        do
        {
            cout << "negativePart: ";
            cin >> inputstr;
            if(inputstr == "h")
            {
                cout << endl <<"Type y for retaining negative part in solutions,, "
                               "otherwise type n." ;
            }
            cout << endl;
        }while(inputstr != "y" && inputstr != "n");

        if(inputstr == "y")
            negativePart = true;
        else
            negativePart = false;
    }

    do
    {
        cout << "Do you supply any extra constant parameter(s)?: " ;
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "Type y for yes.\n"
                            "Type n for no.\n"
                            "GiNaCDE library solve an overdetermined system of algebraic equations only for\n"
                            "constant parameters internally generated. If you wish to solve for any constant\n"
                            "parameters appeared in differential equation type y.";
        }
        cout << endl;
    }while(inputstr != "y" && inputstr != "n");
    if(inputstr == "y")
    {
        do
        {
            cout << "provide extra parameters(s) in curly bracket: ";
            cin >> inputstr;
            if(inputstr == "h" || !is_a<lst>(reader(inputstr)))
            {
                cout << endl << "Provide extra parameters(s) in curly bracket like this : {para1,para2}.";
            }
            cout << endl;
        }while(!is_a<lst>(reader(inputstr)));
        if(nops(ex_to<lst>(reader(inputstr))) != 0)
            paraInDiffSolve = ex_to<lst>(reader(inputstr));
    }

    do
    {
        cout << "Name of output file: ";
        cin >> inputstr;
        if(inputstr == "h")
        {
            cout << endl << "results are saved in this file." ;
        }
        cout << endl;
    }while(inputstr == "h");

    filename = inputstr;


    int ret;

    ret = desolve(inputdiff, {dpndt_var}, method);
    cout<< " "<<endl;
    cout << "Successfully solved." << endl<< endl;

    cin.get(); // consume one newline in buffer
    cin.get();
    return ret;

}
