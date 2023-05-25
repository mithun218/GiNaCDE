
/** @file F_expns_methd.cpp
 *
 *  Implementation of F-expansion and modified F-expansion methods */



#include <ginac/ginac.h>
#include<cln/exception.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "inifcns.h"
#include "simplify.h"
#include "utility.h"
#include "derivative.h"
#include "solve.h"
#include "F_expns_methd.h"
#include "desolve.h"
#include "outform.h"


using namespace std;
using namespace GiNaC;


const symbol _n=symbol("_n");
const ex Fd_=reader("Fd_"), F_ = reader("F_"), F = reader("F");


bool ASolve = false;
bool positivePart = true, negativePart = true;

#ifdef GiNaCDE_gui
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#endif // GiNaCDE_gui

int F_expans::operator()(const ex diffeq, const ex dpndt_varChng, const ex dpndt_var, const ex indpndt_var, const lst twc,
                         const ex tw_coordi, const lst phase, const ex tw_coordiPhase, lst variables, stringstream& solutions,
                         const ex& Nvalue, const int method, const lst& diffDenomSolu, const ex& remainingDiffpart)
{

    string out = outstr("*", 100);
    string out1 = outstr("-", 100);

    //Creating Nvalue symbols
    string symname;
    exvector asym, Asym;
    ex solu_form = _ex0, F_odesubs, odedeg;
    //ex F = reader("F");
    depend(F, {indpndt_var});
    const int NvalueNumer = ex_to<numeric>(numer(Nvalue)).to_int();
    const ex NvalueDenom = ex_to<numeric>(denom(Nvalue));


    if(method == F_expansion || method == mF_expansion)
    {
        if(positivePart)
        {
            if(positivePart)

            for(int i = 0; i < (NvalueNumer) + 1; i++)
            {
                symname = "a_" + to_string(i);
                solu_form = solu_form + pow(F, i/NvalueDenom)*reader(symname);
                asym.push_back(reader(symname));
            }
        }
        if(negativePart)
        {
            for(int i = 1; i < NvalueNumer + 1; i++)
            {
                symname = "b_" + to_string(i);
                solu_form = solu_form + pow(F, _ex_1*i/NvalueDenom)*reader(symname);
                asym.push_back(reader(symname));
            }
        }
    }

    // reset extern variables
    positivePart = true;
    negativePart = true;

    find_ode_order find_ode_order;
    find_ode_order(diffeq);
    int order = *max_element(find_ode_order.order_clt.begin(), find_ode_order.order_clt.end());
    /* relating differentiation of F with F */
    ex F_ode = _ex0;
    exmap diff_F;


    if( method == mF_expansion )
    {

        odedeg = (degAcoeff[0]);
        for(int i = 0; i <= ex_to<numeric>(odedeg).to_int(); i++)
        {
               F_ode = F_ode + pow(F, i)*degAcoeff[i+1];
               Asym.push_back(degAcoeff[i+1]);
        }

        diff_F[Diff(F, indpndt_var, 1)] = F_ode;
    }
    else if( method == F_expansion)
    {
        odedeg = (degAcoeff[0]);
        for(size_t i = 0; i <= ex_to<numeric>(odedeg); i++)
        {
               F_ode = F_ode + pow(F, i)*degAcoeff[i+1];
               Asym.push_back(degAcoeff[i+1]);
        }
        F_ode = sqrt(F_ode);
        diff_F[Diff(F, indpndt_var, 1)] = (F_ode);
    }

    #ifdef GiNaCDE_gui
    gtk_widget_set_sensitive(window,FALSE);
    #endif // GiNaCDE_gui

    const ex C1_=reader("C1_"), C_=reader("C_");

    solutions << dpndt_varChng << " = " << solu_form << ";" << endl;
    cout << dpndt_varChng << " = " << solu_form << ";" << endl;


    solutions << "The first-order nonlinear ODE: ";
    if(output == mathematica)
        solutions <<  "D[F[" <<indpndt_var << "]," << indpndt_var << "] = "<<F_ode <<";" << endl;
    else if(output == maple)
        solutions << "diff(F(" <<indpndt_var << ")," << indpndt_var << ") = "<<F_ode <<";" << endl;
    else
        solutions << "diff(F," << indpndt_var << ",1) = "<<F_ode <<";" << endl;

    solutions << "\n" << out << endl << endl;

    if( order > 1)
    {
        for(int i = 2; i <= order+1; i++)
        {
            diff_F[Diff(F, indpndt_var, i)] = Simplify(expand((pdiff(F_ode, indpndt_var, i-1)).subs(diff_F, subs_options::algebraic)));
        }
    }



    ex solu_form_subs = Simplify2(numer(Simplify(expand(Simplify(expand(evaluate(diffeq.subs(dpndt_varChng == solu_form)))).subs(diff_F, subs_options::algebraic)))));
    lst coeffs;
    //const ex Fd_=reader("Fd_"), F_ = reader("F_");

    solutions << "The system of algebraic equations are: " << endl;
    if(method == F_expansion)
    {
        solu_form_subs=Simplify2(Simplify(expand(solu_form_subs.subs((F_ode.op(0))==pow(Fd_,_ex2),subs_options::algebraic))));

        if(denom(Nvalue) != _ex1)
        {
            solu_form_subs=Simplify2(Simplify(expand(solu_form_subs.subs(F==pow(F_,denom(Nvalue)),subs_options::algebraic))));
            solu_form_subs=solu_form_subs.subs(F_==F,subs_options::algebraic);
        }

        const int Fd_Degree = degree(solu_form_subs,Fd_);

        for(int i = 0; i < degree(solu_form_subs, F) + 1; i++)
        {
            for(int j = 0; j < Fd_Degree + 1; j++)
            {
                coeffs.append((solu_form_subs.coeff(F, i)).coeff(Fd_, j));

                if(denom(Nvalue) != _ex1)
                    solutions <<Diff(F, indpndt_var, 1) << "^" << j <<"*"<< "F^(" << i<<"/"<<denom(Nvalue) << "): "  << (solu_form_subs.coeff(F, i)).coeff(Fd_, j)<< " = 0;" << endl;
                else
                    solutions <<Diff(F, indpndt_var, 1) << "^" << j <<"*"<< F << "^" << i << ": " << (solu_form_subs.coeff(F, i)).coeff(Fd_, j)<< " = 0;" << endl;
            }
        }

    }
    else if( method == mF_expansion )
    {
        if(denom(Nvalue) != _ex1)
        {
            solu_form_subs=Simplify2(expand(solu_form_subs.subs(F==pow(F_,denom(Nvalue)),subs_options::algebraic)));
            solu_form_subs=solu_form_subs.subs(F_==F,subs_options::algebraic);
        }

        for(int i = 0; i < degree(solu_form_subs, F) + 1; i++)
        {
            coeffs.append(solu_form_subs.coeff(F, i));
            if(denom(Nvalue) != _ex1)
                solutions <<"F^(" << i<<"/"<<denom(Nvalue) << "): " << solu_form_subs.coeff(F, i)<< " = 0;" << endl;
            else
                solutions <<"F^" << i << ": " << solu_form_subs.coeff(F, i)<< " = 0;" << endl;
        }

    }

    solutions << " In the following results C_ is an arbitrary constant." << endl;
    solutions << "\n" << out <<  "\n" << endl;


    for(unsigned i = 0; i < asym.size(); i++)
        variables.append(asym[i]);


    ///Assigning traveling wave coordinates
    ex tra_wave_coord;
    if(nops(twc)!=0)
        tra_wave_coord = tw_coordi;
    else
        tra_wave_coord = indpndt_var;

    ex solu_form_F_subs;
    set<lst, ex_is_less>  solu_set_clt;

    lst temvar;
    for( auto it = variables.begin(); it != variables.end(); it++ )
    {
        if( !temvar.has( *it ) )
            temvar.append( *it );
    }
    variables = temvar;
    if( ASolve )
    {
        if(!Asym.empty())
        {
            for(auto it = Asym.begin(); it != Asym.end(); it++)
            {
                exset tem = symbols( *it );
                for( auto it1 = tem.begin(); it1 != tem.end(); it1++ )
                {
                   variables.append(*it1);
                }
            }
        }
    }
    ASolve = false;

    solutions << "solving above system of equations for variables " << variables << "->" << endl << endl;
    cout << "System of equations are solved for the variables " << variables << "....."<< endl;

    solu_set_clt = solve(coeffs, variables );

    if(solu_set_clt.empty())
    {
        solutions << "No solution of the above equations exist;" << endl;
        writetofile(solutions, dpndt_var);

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "No solution exist;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }

    int odetype;
    ex tra_wave_coordSubs,tw_coordiPhaseSubs = _ex0;

    int solNum = 1; // counts solution numbers

    for(set<lst, ex_is_less>::const_iterator it = solu_set_clt.begin(); it != solu_set_clt.end(); it++)
    {
        // solutions of nonlinear algebraic equations //
        solutions << *it << endl;

        replaceI replaceI;
        // In that case of complex diff equation all the parameters are assumed real.
        // Check real and imaginary part have not changed due to present of
        // I in solving parameters.  such as: I*Diff(u,t,1)+p*Diff(u,x,2)+q*u*u*conjugate(u)
        if((remainingDiffpart) != _ex0 &&
            (replaceI((*it)).has(symb_)))
        {
            solutions << "no solutions due to complex parameter;" << endl ;
        }
        else
        {
            try
            {

                F_odesubs = F_ode.subs( *it );

                const ex temF_odesubs = simplify(F_odesubs);
                if(is_polynomial(temF_odesubs,F))
                    F_odesubs = collect(expand(temF_odesubs),F);

                F_odesubs = Simplify2(F_odesubs);

                tra_wave_coordSubs = ((tra_wave_coord.subs(*it)));
                if(nops(phase)!=0)
                    tw_coordiPhaseSubs = tw_coordiPhase.subs(*it);

                solu_form_subs = (solu_form).subs(*it);
                if(solu_form_subs.has(F))
                {
                    solutions << dpndt_varChng << " = "  << solu_form_subs << "," << endl;

                    solutions << "where F is the solution of" << endl;
                    if(output == mathematica)
                        solutions <<  "D[F[" <<indpndt_var << "]," << indpndt_var << "] = "<<F_odesubs <<";" << endl;
                    else if(output == maple)
                        solutions << "diff(F(" <<indpndt_var << ")," << indpndt_var << ") = "<<F_odesubs <<";" << endl;
                    else
                        solutions << "diff(F," << indpndt_var << ",1) = "<<F_odesubs <<";" << endl;;

                    odetype = odeType_check(F_odesubs, F);

                    degAcoeff.remove_all();
                    if(odetype == fracbernouli) // handling y' = A*y^(n/2), y' = sqrt(A*y)
                    {
                        F_odesubs = simplify(F_odesubs);
                        const ex temSubs = Simplify2(expand(F_odesubs.subs(F==pow(Fd_,_ex2))));
                        degAcoeff.append(degree(temSubs, Fd_)/_ex2);
                        degAcoeff.append(Simplify(expand(coeff(temSubs, Fd_,degree(temSubs, Fd_)))));
                    }
                    else if(method == mF_expansion)
                    {
                        degAcoeff.append(degree(F_odesubs, F));
                        for(int i = 0; i <degAcoeff[0]+1; i++)
                            degAcoeff.append(Simplify(expand(coeff(F_odesubs, F,i))));
                    }
                    else
                    {
                        if (!is_a<power>(F_odesubs))
                        {
                            if (is_a<numeric>(F_odesubs))
                            {
                                degAcoeff.append(_ex0);
                                degAcoeff.append(F_odesubs);
                            }
                            else
                            {
                                degAcoeff.append(degree(F_odesubs, F));
                                for (int i = 0; i < degAcoeff[0] + 1; i++)
                                    degAcoeff.append(Simplify(expand(coeff(F_odesubs, F, i))));
                            }
                        }
                        else
                        {
                            degAcoeff.append(degree(F_odesubs.op(0), F));
                            for (int i = 0; i < degAcoeff[0] + 1; i++)
                                degAcoeff.append(Simplify(expand(coeff(F_odesubs.op(0), F, i))));
                        }
                    }

                    lst odeSolu;

                    // getting solutions of first order NLODE //
                    odeSolu = firstOrderDiff_solu(tra_wave_coordSubs, odetype);

                    if(nops(odeSolu) != 0)
                    {
                        solutions << "solution(s) of input Diff. Equ. is (are)=>"<<endl;

                        solutionClt.push_back(lst{});
                        solutionClt[solutionClt.size()-1].append(*it);

                        ex temit;
                        for(auto it1 = odeSolu.begin(); it1 != odeSolu.end(); it1++)
                        {
                            if(is_a<lst>(*it1)) // for handling solutions with conditions
                                temit = (*it1).op(0);
                            else
                                temit = *it1;

                            if(solu_form_subs.has(Diff(wild(0),wild(1),wild(2)))) // handling Diff in solu_form for eF_expansion method
                            {
                                solu_form_F_subs = evaluate(solu_form_subs.subs(F == temit));
                                solu_form_F_subs = Simplify(expand(solu_form_F_subs.subs(xi == tra_wave_coordSubs)));
                            }
                            else if(solu_form_subs.has(F))
                                solu_form_F_subs = solu_form_subs.subs(F == temit);
                            else
                                solu_form_F_subs = solu_form_subs;

                            solu_form_F_subs = (solu_form_F_subs*exp( I*tw_coordiPhaseSubs ));



                            if(nops(diffDenomSolu)!=0)
                            {
                                bool isUndefine = false;
                                for(auto it2 = diffDenomSolu.begin(); it2 != diffDenomSolu.end();)
                                {
                                    if(*it2 == solu_form_F_subs)
                                    {
                                        solutions << "solutions are undefined;"<< endl ;
                                        it2 = diffDenomSolu.end();
                                        isUndefine = true;
                                    }
                                    else
                                    {
                                        it2++;
                                    }
                                }

                                if(!isUndefine)
                                {
                                    if(!is_a<lst>(*it1))
                                    {
                                        solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << solu_form_F_subs <<";"<< endl;
                                        solNum = solNum + 1;

                                        solutionClt[solutionClt.size()-1].append(dpndt_var  ==  solu_form_F_subs);
                                    }
                                    else // for handling solutions with conditions
                                    {
                                        lst temList = {};
                                        solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << solu_form_F_subs <<" (with condition(s) ";
                                        temList.append(dpndt_var  ==  solu_form_F_subs);

                                        for(size_t i=1; i<nops(*it1); i++)
                                        {
                                            solutions<<(*it1).op(i)<<", ";
                                            temList.append((*it1).op(i));
                                        }
                                        solutions<<");"<<endl;

                                        solNum = solNum + 1;

                                        solutionClt[solutionClt.size()-1].append(temList);
                                    }
                                }

                            }
                            else
                            {
                                if(!is_a<lst>(*it1))
                                {
                                    solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << solu_form_F_subs <<";"<< endl;
                                    solNum = solNum + 1;

                                    solutionClt[solutionClt.size()-1].append(dpndt_var  ==  solu_form_F_subs);
                                }
                                else
                                {
                                    solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << solu_form_F_subs <<" (with condition(s) ";

                                    lst temList = {};
                                    temList.append(dpndt_var  ==  solu_form_F_subs);

                                    for(size_t i=1; i<nops(*it1); i++)
                                    {
                                        solutions<<(*it1).op(i)<<", ";
                                        temList.append((*it1).op(i));
                                    }
                                    solutions<<");"<<endl;
                                    solNum = solNum + 1;

                                    solutionClt[solutionClt.size()-1].append(temList);

                                }

                            }


                        }

                    }


                }
                else
                {
                    solutions << "solution(s) of input Diff. Equ. is (are)=>"<<endl;

                    if(nops(diffDenomSolu)!=0)
                    {
                        bool isUndefine = false;
                        for(auto it1 = diffDenomSolu.begin(); it1 != diffDenomSolu.end();)
                        {
                            if(*it1 == solu_form_subs)
                            {
                                solutions << "solutions are undefined;"<< endl ;
                                it1 = diffDenomSolu.end();
                                isUndefine = true;
                            }
                            else
                            {
                                it1++;
                            }
                        }
                        if(!isUndefine)
                        {
                            solutions <<"solution #"<<solNum<<"  "<< dpndt_var << " = "  << solu_form_subs*exp( I*tw_coordiPhaseSubs )<<";"<< endl;

                            solutionClt.push_back(lst{});
                            solutionClt[solutionClt.size()-1].append(*it);
                            solutionClt[solutionClt.size()-1].append(dpndt_var  ==  solu_form_subs*exp( I*tw_coordiPhaseSubs));
                        }

                    }
                    else
                    {
                        solutions <<"solution #"<<solNum<<"  "<< dpndt_var << " = "  << solu_form_subs*exp( I*tw_coordiPhaseSubs ) <<";"<< endl;
                        solNum = solNum + 1;
                        solutionClt.push_back(lst{});
                        solutionClt[solutionClt.size()-1].append(*it);
                        solutionClt[solutionClt.size()-1].append(dpndt_var  ==  solu_form_subs*exp( I*tw_coordiPhaseSubs));
                    }

                }

            }
            catch(GiNaC::pole_error){solutions<<"GiNaC::pole_error" <<endl;}
            catch (cln::runtime_exception){solutions<<"cln::runtime_exception" <<endl;}

            catch(std::invalid_argument){solutions<<"std::invalid_argument" <<endl;}
            catch(std::out_of_range){solutions<<"std::out_of_range" <<endl;}
            //catch(std::runtime_error){solutions<<"std::runtime_error" <<endl;}
            //catch(std::range_error){solutions<<"std::range_error" <<endl;}
            //catch(std::logic_error){solutions<<"std::logic_error" <<endl;}
           // catch(std::domain_error){solutions<<"std::domain_error" <<endl;}


        }

        solutions<< endl;

    }
    solutions <<  "\n" << out << endl;
    solutions << out << "\n" << endl;



    const std::chrono::time_point<std::chrono::system_clock> endTime = chrono::high_resolution_clock::now();
    auto dur = endTime-beginTime;
    cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;
    solutions << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;
    writetofile(solutions, dpndt_var);

    #ifdef GiNaCDE_gui
    stringstream temstr;
    temstr << "The results are written in " << filename <<" file.             " << "Time:"
    << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" ;

    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, &temstr.str()[0]);
    //gtk_show_uri_on_window(GTK_WINDOW(window),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
    
    gtk_widget_set_sensitive(window,TRUE);
    #endif // GiNaCDE_gui

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////



