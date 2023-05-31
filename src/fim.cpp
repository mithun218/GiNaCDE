
/** @file fim.cpp
 *
 *  Implementation of first integral method*/



#include <ginac/ginac.h>
#include<cln/exception.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "inifcns.h"
#include "simplify.h"
#include "derivative.h"
#include "desolve.h"
#include "utility.h"
#include "solve.h"
#include "integrate.h"
#include "fim.h"
#include "outform.h"

using namespace std;
using namespace GiNaC;


const ex Y_=reader("Y_"), X_=reader("X_");

/////////////////////////////////////////////////////////////
ex degreePartiVar::operator()(const ex & e)
{
    if( is_a<GiNaC::function>(e) )
        return e;
    else if(is_a<power>(e) && e.op(0) == partiVar)
    {
        degreeClt.push_back(e.op(1));
        return e;
    }
    else if( e == partiVar )
    {
        degreeClt.push_back(_ex1);
        return e;
    }
    else
    {
        return e.map(*this);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
ex fim::removeex(const ex& _expr1, const ex& _expr2)
{
    ex removed = _ex0;
    if(is_a<add>(_expr1))
    {
        for(unsigned i = 0; i < nops(_expr1); i++)
        {
            if(!_expr1.op(i).has(_expr2))
            {
                removed = removed + _expr1.op(i);
            }
        }
    }

    return _ex_1*removed;
}
/////////////////////////////////////////////////////////////////////

/** Solutions of linear and quadratic polynomial equation used in fim **/
lst quadraticFim(const ex & equ_, const ex& var_)
{
    lst quasolu;

    if(equ_.degree(var_) == 0)
        return quasolu;

    if(equ_.degree(var_) == 1)
        quasolu.append(var_ == simplify(Simplify(_ex_1*equ_.coeff(var_, 0)))/(Simplify(equ_.coeff(var_, 1))));
    else
    {
        ex _a, _b, _c;

        _a = Simplify(equ_.coeff(var_, 2));
        _b = Simplify(equ_.coeff(var_, 1));
        _c = Simplify(equ_.coeff(var_, 0));
        ex tem = Simplify(expand(_b*_b - 4*_a*_c));
        if(tem == _ex0)
        {
            quasolu.append(var_ == Simplify((_ex1_2*(-_b))/_a));
        }
        else
        {
            quasolu.append(var_ == (Simplify((-_b/(2*_a) + sqrt(Simplify(expand(tem/(4*pow(_a,2)))))))));
            quasolu.append(var_ == (Simplify((-_b/(2*_a) - sqrt(Simplify(expand(tem/(4*pow(_a,2)))))))));
        }
    }

    for(unsigned i = 0; i < nops(quasolu); i++)
    {
        if(quasolu.op(i).rhs().info(info_flags::positive) && quasolu.op(i).rhs() < Gtolerance)
            quasolu[i] = {var_ == _ex0};
        else if(quasolu.op(i).rhs().info(info_flags::negative) && -quasolu.op(i).rhs() < Gtolerance)
            quasolu[i] = {var_ == _ex0};
    }

    return quasolu;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
int fim::operator()(const ex diffeq, const ex dpndt_varChng, const ex dpndt_var, const ex indpndt_var, const lst twc,
                    const ex tw_coordi, const lst phase, const ex tw_coordiPhase, lst variables, stringstream& solutions,
                    const lst& diffDenomSolu, const ex& remainingDiffpart)
{
    #ifdef GiNaCDE_gui
    gtk_widget_set_sensitive(window,FALSE);
    #endif // GiNaCDE_gui

    string out = outstr("*", 100);

    int Nvalue;
    if(NValue <= 0)
        Nvalue = 1;
    else
        {Nvalue = ex_to<numeric>(NValue).to_int();NValue=_ex0;}

    NValue = 0;

    solutions << "The value of N is: " << Nvalue << ";" << endl;
    solutions << "\n" << out << endl << endl;

    const int ordr = order(diffeq);
    if(ordr > 2)
    {
        solutions << "Evaluation stop: order of ode should be 2;" << endl;
        cout << "Evaluation stop: order of ode should be 2;" << endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: order of ode should be 2;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }

    exmap Null;
    const lst hdifftrml = collectAllCoeff(diffeq, {Diff(dpndt_varChng,indpndt_var,ordr)},false,Null);

    if(nops(hdifftrml)>1 || nops(hdifftrml)==0)
    {
       solutions << "Evaluation stop: unsupported ode;" << endl;
       cout << "Evaluation stop: unsupported ode;" <<hdifftrml<<"  "<<diffeq<< endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported ode;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }

    ex hdifftrm = hdifftrml[0];
    ex lhsodetrm;

    lhsodetrm = Simplify(expand(removeex(diffeq, Diff(dpndt_varChng,indpndt_var,ordr))));

    ///Assigning traveling wave coordinates
    ex tra_wave_coord;
    if(nops(twc)!=0)
        tra_wave_coord = tw_coordi;
    else
        tra_wave_coord = indpndt_var;


    //const ex Y_=reader("Y_"), X_=reader("X_");
    const ex g_=reader("g_"), h_=reader("h_");

    depend(Y_, {X_});
    depend(g_, {Y_});
    depend(h_, {Y_});

    lhsodetrm = lhsodetrm.subs(lst{Diff(wild(0), wild(1), wild(2)) == Y_});
    lhsodetrm = lhsodetrm.subs(lst{dpndt_varChng == X_});
    hdifftrm = hdifftrm.subs(lst{dpndt_varChng == X_});

    if(lhsodetrm.degree(Y_) - 1 > Nvalue)
    {
        solutions << "Evaluation stop: provide minimum Value of N: " << lhsodetrm.degree(Y_) - 1 << endl;
        cout << "Evaluation stop: provide minimum Value of N: " << lhsodetrm.degree(Y_) - 1 << endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "No solution exist;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }

    if(lhsodetrm.has(indpndt_var))
    {
        solutions << "Evaluation stop: unsupported ode;" << endl;
        cout << "Evaluation stop: unsupported ode;" << endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported ode;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }


    ex rhstrm = _ex0, lhstrm = _ex0;
    string symname;
    exvector asym, arbtCnst;

    ex exWtoutX = _ex1, exWtX = _ex1;;

    if(lhsodetrm.degree(Y_)==_ex2 && hdifftrm.has(X_))
    {
        hdifftrm = Simplify(collect_common_factors(Factor(hdifftrm)));


        if(is_a<mul>(hdifftrm))
        {
            for(size_t i = 0; i<nops(hdifftrm); i++)
            {
                if(hdifftrm.op(i).has(X_))
                {
                    exWtX = exWtX*hdifftrm.op(i);
                }
                else
                {
                    exWtoutX = exWtoutX*hdifftrm.op(i);
                }
            }
            hdifftrm = exWtX;
        }

        solutions << "We make the transformation, d xi = (" <<hdifftrm<<")*d eta to avoid singularity "<<hdifftrm <<" = 0 temporarily."<< endl;
        cout <<"Now, we make transformation, d xi = (" <<hdifftrm<<")*d eta, to avoid singularity "<<hdifftrm <<" = 0, temporarily." << endl;
        lhsodetrm = collect(Simplify(expand(lhsodetrm/exWtoutX)),Y_);

        if(output==maple)
        {
            solutions << "Let "<<dpndt_varChng<<" = X_, diff("<<dpndt_varChng<<",eta) = Y_*("<<hdifftrm<<"), then we get\n diff(X_,eta) = Y_*("<<hdifftrm<<"),\n"<<
                         "diff(Y_,eta) = "<<lhsodetrm <<","<< endl << endl;
        }
        else if(output==mathematica)
        {
            solutions << "Let "<<dpndt_varChng<<" = X_, D["<<dpndt_varChng<<",eta] = Y_*("<<hdifftrm<<"), then we get\n D[X_,eta] = Y_*("<<hdifftrm<<"),\n"<<
                         "D[Y_,eta] = "<<lhsodetrm <<","<< endl << endl;
        }
        else
        {
            solutions << "Let "<<dpndt_varChng<<" = X_, Diff("<<dpndt_varChng<<",eta, 1) = Y_*("<<hdifftrm<<"), then we get\n Diff(X_,eta, 1) = Y_*("<<hdifftrm<<"),\n"<<
                         "Diff(Y_,eta, 1) = "<<lhsodetrm <<","<< endl << endl;
        }

        lhsodetrm = Simplify(expand(lhsodetrm));
    }
    else
        lhsodetrm = Simplify(expand(lhsodetrm/hdifftrm));

    for(int i = 0; i < Nvalue+1; i++)
    {
        symname = "a_" + to_string(i);
        asym.push_back(reader(symname));
        depend(asym[i], {X_});
        lhstrm = lhstrm + Diff(asym[i], X_, 1)*pow(Y_, i+1)*exWtX ;
        rhstrm = rhstrm + (g_ + h_*Y_)*pow(Y_, i)*asym[i]- i*pow(Y_, i-1)*asym[i]*(lhsodetrm);
    }
    lhstrm = lhstrm + asym[1]*lhsodetrm.coeff(Y_,0);
    rhstrm = rhstrm + asym[1]*lhsodetrm.coeff(Y_,0);

    const ex denoma  = denom(Simplify(expand(simplify(expand(lhstrm - rhstrm)))));
    ex finalequ = Simplify(expand(simplify(expand((lhstrm - rhstrm)*denoma))));
    lhstrm = Simplify(expand(simplify(expand((lhstrm)*denoma))));
    rhstrm = Simplify(expand(simplify(expand((rhstrm)*denoma))));
    lst substit;

    const int yddegree = finalequ.degree(Y_);

    substit.append(asym[Nvalue] == _ex1);

    lst asymlst={};
    ex lhstermTem;
    for(auto it=asym.begin(); it!=asym.end();it++)
        asymlst.append(*it);
    solutions<<diffformchange(lhstrm,asymlst,{X_})<<" = "<<rhstrm<<endl;

    solutions << "Comparing the coefficients of " << Y_<<"^i (i =" << yddegree << " .., 0) in both sides, we have" << endl << endl;

    for(int i = yddegree; i >= 0; i--)
    {
        solutions <<diffformchange(Simplify(lhstrm.coeff(Y_, i)), asymlst,{X_}) << " = " << Simplify(rhstrm).coeff(Y_, i) <<","<< endl;
    }


    ex ydegCoeff = (finalequ.coeff(Y_, yddegree));
    ydegCoeff = Simplify(expand(evaluate(ydegCoeff.subs(asym[Nvalue]==_ex1))));

    ex hValue = ydegCoeff+h_;
    if(hValue != _ex0 && ydegCoeff.has(h_))
    {
        const exsetlst solutemp = solve({ydegCoeff},{h_});

        if(!solutemp.empty())
            hValue = (*solutemp.begin())[0].rhs();

    }

    if(hValue != _ex0 && !is_polynomial(hValue,X_))
    {
        solutions << "Evaluation stop: unsupported ode;" << endl;
        cout << "Evaluation stop: unsupported ode;" << endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported ode;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui

        return -1;
    }

    solutions << "assuming " << asym[Nvalue]<<" = 1, in first equation, we get h_ = "<<hValue <<";" << endl;
    substit.append(h_ == hValue);

    /// Balancing degree of X_ in the above system of equations
    degreePartiVar degreePartiVar(X_);
    bool isBalance = true, isNegative=false;
    int hdegree = 6;
    vector<vector<int>> balancedDegree;

    if(Nvalue == 1)
    {
        std::vector<int> lhsSubs, rhsSubs;
        ex a0Deg = reader("a0Deg"),gDeg = reader("gDeg");
        ex lhstrmSubs = lhstrm.subs(substit); ex rhstrmSubs = rhstrm.subs(substit);
        lhstrmSubs = Simplify(expand(evaluate(lhstrmSubs.subs(lst{asym[0]==pow(X_,a0Deg),g_==pow(X_,gDeg)}))));
        rhstrmSubs = Simplify(expand(evaluate(rhstrmSubs.subs(lst{asym[0]==pow(X_,a0Deg),g_==pow(X_,gDeg)}))));

        exvector tem;
        ex temex;
        vector<vector<ex>> lhsDegClt, rhsDegClt;
        for(int i = yddegree-1; i >= 0; i--)
        {
            degreePartiVar.degreeClt.clear();
            degreePartiVar(Simplify(expand(lhstrmSubs.coeff(Y_, i))));

            temex = lhstrmSubs.coeff(Y_, i);
            tem = degreePartiVar.degreeClt;
            if(is_a<add>(temex))
            {
                for(size_t i = 0; i<nops(temex);)
                {
                    if(!temex.has(X_))
                    {
                        tem.push_back(_ex0);
                        i = nops(temex);
                    }
                    else
                    {
                        i++;
                    }
                }
                lhsDegClt.push_back(tem);
            }
            else
            {
                if(!temex.has(X_))
                {
                    tem.push_back(_ex0);
                }
                lhsDegClt.push_back(tem);
            }


            degreePartiVar.degreeClt.clear();
            degreePartiVar(Simplify(expand(rhstrmSubs.coeff(Y_, i))));

            temex = rhstrmSubs.coeff(Y_, i);
            tem = degreePartiVar.degreeClt;
            if(is_a<add>(temex))
            {
                for(size_t i = 0; i<nops(temex);)
                {
                    if(!temex.has(X_))
                    {
                        tem.push_back(_ex0);
                        i = nops(temex);
                    }
                    else
                    {
                        i++;
                    }
                }
                rhsDegClt.push_back(tem);
            }
            else
            {
                if(!temex.has(X_))
                {
                    tem.push_back(_ex0);
                }
                rhsDegClt.push_back(tem);
            }

        }

        for(int a0Degree = hdegree;a0Degree>=0;a0Degree--)
        {
            for(int gDegree = hdegree;gDegree>=0;gDegree--)
            {
                isBalance = true;
                isNegative = false;
                for(unsigned i = 0; i < lhsDegClt.size(); i++)
                {
                    lhsSubs.clear();
                    rhsSubs.clear();
                    for( auto it = lhsDegClt[i].begin(); it != lhsDegClt[i].end(); it++ )
                    {
                        if(ex_to<numeric>((*it).subs(lst{a0Deg==a0Degree,gDeg==gDegree})).to_int()<0)
                            isNegative=true;
                        else
                            lhsSubs.push_back(ex_to<numeric>((*it).subs(lst{a0Deg==a0Degree,gDeg==gDegree})).to_int());
                    }

                    for( auto it = rhsDegClt[i].begin(); it != rhsDegClt[i].end(); it++ )
                    {
                        if(ex_to<numeric>((*it).subs(lst{a0Deg==a0Degree,gDeg==gDegree})).to_int()<0)
                            isNegative=true;
                        else
                            rhsSubs.push_back(ex_to<numeric>((*it).subs(lst{a0Deg==a0Degree,gDeg==gDegree})).to_int());
                    }

                    if(*max_element(lhsSubs.begin(),lhsSubs.end())!=(*max_element(rhsSubs.begin(),rhsSubs.end())) || isNegative)
                    {
                        isBalance = false;
                    }

                }

                if(isBalance && !isNegative)
                    balancedDegree.push_back({a0Degree,gDegree});

            }

        }
    }

    else if(Nvalue == 2)
    {
        exvector lhstemvector, rhstemvector;
        std::vector<int> lhsSubs, rhsSubs;
        ex a1Deg = reader("a1Deg"),a0Deg = reader("a0Deg"),gDeg = reader("gDeg"),
        YlhsEq, YrhsEq, Y2EqXDegree, Y1EqXDegree,Y0EqXDegree;
        ex lhstrmSubs = lhstrm.subs(substit); ex rhstrmSubs = rhstrm.subs(substit);
        lhstrmSubs = Simplify(expand(evaluate(lhstrmSubs.subs(lst{asym[1]==pow(X_,a1Deg),asym[0]==pow(X_,a0Deg),g_==pow(X_,gDeg)}))));
        rhstrmSubs = Simplify(expand(evaluate(rhstrmSubs.subs(lst{asym[1]==pow(X_,a1Deg),asym[0]==pow(X_,a0Deg),g_==pow(X_,gDeg)}))));

        exvector tem;
        ex temex;
        vector<vector<ex>> lhsDegClt, rhsDegClt;
        for(int i = yddegree-1; i >= 0; i--)
        {
            degreePartiVar.degreeClt.clear();
            degreePartiVar(Simplify(expand(lhstrmSubs.coeff(Y_, i))));

            temex = lhstrmSubs.coeff(Y_, i);
            tem = degreePartiVar.degreeClt;
            if(is_a<add>(temex))
            {
                for(size_t i = 0; i<nops(temex);)
                {
                    if(!temex.has(X_))
                    {
                        tem.push_back(_ex0);
                        i = nops(temex);
                    }
                    else
                    {
                        i++;
                    }
                }
                lhsDegClt.push_back(tem);
            }
            else
            {
                if(!temex.has(X_))
                {
                    tem.push_back(_ex0);
                }
                lhsDegClt.push_back(tem);
            }


            degreePartiVar.degreeClt.clear();
            degreePartiVar(Simplify(expand(rhstrmSubs.coeff(Y_, i))));

            temex = rhstrmSubs.coeff(Y_, i);
            tem = degreePartiVar.degreeClt;
            if(is_a<add>(rhstrmSubs.coeff(Y_, i)))
            {
                for(size_t i = 0; i<nops(temex);)
                {
                    if(!temex.has(X_))
                    {
                        tem.push_back(_ex0);
                        i = nops(temex);
                    }
                    else
                    {
                        i++;
                    }
                }
                rhsDegClt.push_back(tem);
            }
            else
            {
                if(!temex.has(X_))
                {
                    tem.push_back(_ex0);
                }
                rhsDegClt.push_back(tem);
            }

        }


        for(int a1Degree = hdegree;a1Degree>=0;a1Degree--)
        {
            for(int a0Degree = hdegree;a0Degree>=0;a0Degree--)
            {
                for(int gDegree = hdegree;gDegree>=0;gDegree--)
                {
                    isBalance = true;
                    isNegative = false;

                    for(unsigned i = 0; i < lhsDegClt.size(); i++)
                    {
                        lhsSubs.clear();
                        rhsSubs.clear();
                        for( auto it = lhsDegClt[i].begin(); it != lhsDegClt[i].end(); it++ )
                        {

                            if(ex_to<numeric>((*it).subs(lst{a1Deg==a1Degree,a0Deg==a0Degree,gDeg==gDegree})).to_int()<0)
                                isNegative=true;
                            else
                                lhsSubs.push_back(ex_to<numeric>((*it).subs(lst{a1Deg==a1Degree,a0Deg==a0Degree,gDeg==gDegree})).to_int());
                        }

                        for( auto it = rhsDegClt[i].begin(); it != rhsDegClt[i].end(); it++ )
                        {

                            if(ex_to<numeric>((*it).subs(lst{a1Deg==a1Degree,a0Deg==a0Degree,gDeg==gDegree})).to_int()<0)
                                isNegative=true;
                            else
                                rhsSubs.push_back(ex_to<numeric>((*it).subs(lst{a1Deg==a1Degree,a0Deg==a0Degree,gDeg==gDegree})).to_int());
                        }

                        if(*max_element(lhsSubs.begin(),lhsSubs.end())!=(*max_element(rhsSubs.begin(),rhsSubs.end())) || isNegative)
                        {
                            isBalance = false;
                        }

                    }

                    if(isBalance && !isNegative)
                        balancedDegree.push_back({a0Degree,a1Degree,gDegree});

                }

            }

        }
    }


    if(balancedDegree.empty())
    {
        cout << "Evaluation stop: Balancing degree of X_ is failure." << endl;
        solutions << "Evaluation stop: Balancing degree of X_ is failure." << endl;

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Balancing degree of X_ is failure.");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        gtk_widget_set_sensitive(window,TRUE);
        #endif // GiNaCDE_gui
        return -1;
    }
    else
    {
        if(Nvalue==1)
        {
            solutions << "Balancing degrees of X_ we get, degrees of (a_0, g_) = ";
            cout << "Balancing degrees of X_ we get, degrees of (a_0, g_) = ";
        }
        else if(Nvalue==2)
        {
            solutions << "Balancing degrees of X_ we get, degrees of (a_0, a_1, g_) = ";
            cout << "Balancing degrees of X_ we get, degrees of (a_0, a_1, g_) = ";
        }


        for(unsigned i = 0;i<balancedDegree.size();i++)
        {
            solutions<<"(";
            cout<<"(";
            for(unsigned j = 0;j<balancedDegree[i].size();j++)
            {
                solutions<<balancedDegree[i][j];
                cout<<balancedDegree[i][j];
                if(j!=balancedDegree[i].size()-1)
                {
                    solutions<<", ";
                    cout<<", ";
                }
            }
            solutions<<")";
            cout<<"), ";
            if(i!=balancedDegree.size()-1)
                solutions<<", ";
            else
                solutions<<endl<<endl;
        }
        cout<<endl;
    }

    lst sysequ,temvariables;
    ex finalequSubs;
    const string outBkslash = outstr("\\",100);
    int solNum = 1; // counts solution numbers
    for(unsigned blnci1 = 0;blnci1<balancedDegree.size();blnci1++)
    {
        if(Nvalue==1)
            solutions<<"                                //////////Degrees of (a_0, g_) = (";
        else if(Nvalue==2)
            solutions<<"                                //////////Degrees of (a_0, a_1, g_) = (";
        for(unsigned blnci2 = 0;blnci2<balancedDegree[blnci1].size();blnci2++)
        {
            solutions<<balancedDegree[blnci1][blnci2];
            if(blnci2!=balancedDegree[blnci1].size()-1)
                solutions<<", ";
        }
        solutions<<")//////////"<<endl;
        solutions<<"                 "<<outBkslash<<endl;

        substit.remove_all();
        substit.append(h_ == hValue);
        substit.append(asym[Nvalue] == _ex1);
        temvariables.remove_all();
        temvariables=variables;
        if(Nvalue==1)
        {
            ex _a0ex = _ex0;
            for(int i = 0; i <= balancedDegree[blnci1][0]; i++)
            {
                symname = "a_0" + to_string(i);
                temvariables.append(reader(symname));
                _a0ex = _a0ex + pow(X_, i)*reader(symname);
            }
            substit.append(asym[0] == _a0ex);

            ex _gex = _ex0;
            for(int i = 0; i <= balancedDegree[blnci1][1]; i++)
            {
                symname = "g_" + to_string(i);
                temvariables.append(reader(symname));
                _gex = _gex + pow(X_, i)*reader(symname);
            }
            substit.append(g_ == _gex);

            solutions << "hence a_0 = " << _a0ex << "," << endl;
            solutions << "      g_ = " << _gex << ";" << endl;
        }
        else if(Nvalue==2)
        {
            ex _a0ex = _ex0;
            for(int i = 0; i <= balancedDegree[blnci1][0]; i++)
            {
                symname = "a_0" + to_string(i);
                temvariables.append(reader(symname));
                _a0ex = _a0ex + pow(X_, i)*reader(symname);
            }
            substit.append(asym[0] == _a0ex);

            ex _a1ex = _ex0;
            for(int i = 0; i <= balancedDegree[blnci1][1]; i++)
            {
                symname = "a_1" + to_string(i);
                temvariables.append(reader(symname));
                _a1ex = _a1ex + pow(X_, i)*reader(symname);
            }
            substit.append(asym[1] == _a1ex);

            ex _gex = _ex0;
            for(int i = 0; i <= balancedDegree[blnci1][2]; i++)
            {
                symname = "g_" + to_string(i);
                temvariables.append(reader(symname));
                _gex = _gex + pow(X_, i)*reader(symname);
            }
            substit.append(g_ == _gex);

            solutions << "hence a_0 = " << _a0ex << "," << endl;
            solutions << "      a_1 = " << _a1ex << "," << endl;
            solutions << "      g_ = " << _gex << ";" << endl;
        }




        sysequ.remove_all();
        finalequSubs = Simplify(expand(evaluate(finalequ.subs(substit))));

        solutions << "Substituting "<<asym[Nvalue]<<",.. a_0, g_ into all equation and setting all the coefficients of powers of  X_ to zero,"<< endl;
        for(int i = 0; i < finalequSubs.degree(Y_) + 1; i++ )
        {
            for(int j = 0; j < finalequSubs.degree(X_) + 1; j++ )
            {
                sysequ.append((finalequSubs.coeff(Y_, i)).coeff(X_, j));
                solutions <<Y_<<"^"<<i<<"*"<< X_ << "^" << j << ": " << finalequSubs.coeff(Y_, i).coeff(X_, j) << " = 0," << endl;
            }
        }

        solutions << " In the following results C_ is an arbitrary constant." << endl;
        solutions << "\n" << out <<  "\n" << endl;

        lst temvar;
        for( auto it = temvariables.begin(); it != temvariables.end(); it++ )
        {
            if( !temvar.has( *it ) )
                temvar.append( *it );
        }

        temvariables = temvar;

        solutions << "solving above system of equations for variables " << temvariables << "->" << endl << endl;
        cout << "System of equations are solved for the variables " << temvariables << "....."<< endl;
        exsetlst solu = solve(sysequ, temvariables);

        if(solu.empty())
        {
            solutions << "No solution of the system of algebraic equations exist;" << endl;
            cout << "No solution of system of algebraic equations exist;" << endl;
        }
        else
        {
            ///converting 1st order ode
            ex fstordrode = _ex0;
            unsigned po = 0;
            for(exvector::const_iterator it = asym.begin(); it != asym.end(); it++)
            {

                fstordrode = fstordrode + pow(Y_, po)*(*it);
                po = po + 1;
            }

            fstordrode = Simplify(expand(fstordrode.subs(substit)));

            ex solu_form_F_subs1, temfstordrode;
            int odetype;
            ex fstordrodeSubs,tra_wave_coordSubs,tw_coordiPhaseSubs=_ex0;


            for(set<lst, ex_is_less>::const_iterator it = solu.begin(); it != solu.end(); it++)
            {
                solutions << *it << endl;

                replaceI replaceI;
                // In that case of complex diff equation all the parameters are assumed real.
                // Check real and imaginary part have not changed due to present of
                // I in solving parameters.  such as: I*Diff(u,t,1)+p*Diff(u,x,2)+q*u*u*conjugate(u)
                if((remainingDiffpart) != _ex0 &&
                ( replaceI((*it)).has(symb_)))
                {
                    solutions << "no solutions due to complex parameter;" << endl ;
                }
                else
                {
                    try
                    {
                        temfstordrode = fstordrode.subs(*it);
                        string  fstordrodeStr = diffformchange(collect(Simplify(expand(((temfstordrode.subs((lst){X_ == dpndt_varChng}))).subs(Y_ == Diff(dpndt_varChng, indpndt_var, 1)))), Diff(dpndt_varChng, indpndt_var, 1)), {dpndt_varChng}, {indpndt_var});
                        solutions << fstordrodeStr << " = 0," << endl;


                        // first integral forms are collected by this lst variable //
                        lst fstordrodeRhs = {};


                        if(Nvalue>1)
                        {

                            if(temfstordrode.degree(Y_)>1)
                            {

                                // factorizing first integral form //
                                const lst fstordrode_solu = quadraticFim(temfstordrode,Y_);

                                ex factorex = _ex1;

                                if(nops(fstordrode_solu))
                                {
                                    ex fstordrodeSimp, prev;
                                    solutions<< "after factorization of above nonlinear ODE we get: "<<endl;

                                    for(auto it1 = fstordrode_solu.begin(); it1 != fstordrode_solu.end(); it1++)
                                    {

                                        const ex temPolychk = expand(simplify(expand((((*it1).rhs())))));
                                        if(temPolychk.is_polynomial(X_))
                                            fstordrodeRhs.append(temPolychk);
                                        else
                                            fstordrodeRhs.append(Simplify(expand((((*it1).rhs())))));
                                        factorex = factorex*((*it1).lhs()-Simplify(((((*it1).rhs())))));

                                    }

                                    fstordrodeStr = diffformchange(((((factorex.subs(lst{X_ == dpndt_varChng}))).subs(Y_ == Diff(dpndt_varChng, indpndt_var, 1)))), {dpndt_varChng}, {indpndt_var});
                                    solutions << fstordrodeStr << " = 0," << endl;
                                    solutions<<"solutions of each factor in above equation will be determined."<<endl;
                                }
                            }
                            else if(is_polynomial((temfstordrode), X_))
                            {
                                fstordrodeRhs.append(Simplify(expand(collect((temfstordrode- Y_)*_ex_1,X_))));
                            }
                        }
                        else
                            fstordrodeRhs.append(Simplify(expand(collect((temfstordrode - Y_)*_ex_1,X_))));

                        tra_wave_coordSubs = ((tra_wave_coord.subs(*it)));

                        if(nops(phase)!=0)
                            tw_coordiPhaseSubs = tw_coordiPhase.subs(*it);


                        degAcoeff.remove_all();
                        if(nops(fstordrodeRhs) != 0)
                        {
                            for(auto it1 = fstordrodeRhs.begin(); it1 != fstordrodeRhs.end(); it1++ )
                            {
                                const ex temSimp = Simplify2(*it1);

                                odetype = odeType_check(temSimp,X_ );
                                if(is_polynomial(temSimp, X_)) //
                                {
                                    degAcoeff.append(degree(temSimp, X_));
                                    for(unsigned i = 0; i <degAcoeff[0]+1; i++)
                                        degAcoeff.append(coeff(temSimp, X_,i));
                                }
                                else if(odetype == fracbernouli) // handling y' = A*y^(n/2), y' = sqrt(A*y)
                                {
                                    ex temEx = simplify(temSimp),Fd_;
                                    const ex temSubs = Simplify2(Simplify(expand(temEx.subs(X_==pow(Fd_,_ex2)))));
                                    degAcoeff.append(degree(temSubs, Fd_)/_ex2);
                                    degAcoeff.append(Simplify(expand(coeff(temSubs, Fd_,degree(temSubs, Fd_)))));
                                }
                                else if(odetype != general)
                                {
                                    degAcoeff.append(degree((temSimp).op(0), X_));
                                    for (int i = 0; i < degAcoeff[0] + 1; i++)
                                        degAcoeff.append(Simplify(expand(coeff((temSimp).op(0), X_, i))));
                                }


                                // getting solutions of first order NLODE //
                                lst odeSolu = firstOrderDiff_solu(tra_wave_coordSubs, odetype);
                                degAcoeff.remove_all();


                                if(nops(odeSolu) != 0)
                                {
                                    solutions << "solution(s) of input Diff. Equ. is (are)=>"<<endl;

                                    solutionClt.push_back(lst{});
                                    solutionClt[solutionClt.size()-1].append(*it);

                                    ex finaSolu;
                                    for(auto it2 = odeSolu.begin(); it2 != odeSolu.end(); it2++)
                                    {
                                        ex temit;
                                        if(is_a<lst>(*it2)) // for handling solutions with conditions
                                            temit = (*it2).op(0);
                                        else
                                            temit = *it2;

                                        if(nops(diffDenomSolu)!=0)
                                        {
                                            bool isUndefine = false;
                                            for(auto it3 = diffDenomSolu.begin(); it3 != diffDenomSolu.end();)
                                            {
                                                if(*it3 == (temit)*exp( I*tw_coordiPhaseSubs ))
                                                {
                                                    solutions << "solutions are undefined;"<< endl << endl;
                                                    it3 = diffDenomSolu.end();
                                                    isUndefine = true;
                                                }
                                                else
                                                    it3++;
                                            }

                                            if(!isUndefine)
                                            {
                                                if(!is_a<lst>(*it2))
                                                {
                                                    solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << (temit)*exp( I*tw_coordiPhaseSubs ) <<";"<< endl;
                                                    solNum = solNum + 1;

                                                    solutionClt[solutionClt.size()-1].append(dpndt_var  ==  (temit)*exp( I*tw_coordiPhaseSubs ));
                                                }
                                                else // for handling solutions with conditions
                                                {
                                                    solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << (temit)*exp( I*tw_coordiPhaseSubs ) <<" (with condition(s) ";

                                                    lst temList = {};
                                                    temList.append(dpndt_var  ==  (temit)*exp( I*tw_coordiPhaseSubs ));
                                                    for(size_t i=1; i<nops(*it2); i++)
                                                    {
                                                        solutions<<(*it2).op(i)<<", ";
                                                        temList.append((*it).op(i));
                                                    }

                                                    solutions<<");"<<endl;
                                                    solNum = solNum + 1;

                                                    solutionClt[solutionClt.size()-1].append(temList);

                                                }
                                            }
                                        }
                                        else
                                        {
                                            if(!is_a<lst>(*it2))
                                            {
                                                solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << (temit)*exp( I*tw_coordiPhaseSubs ) <<";"<< endl;
                                                solNum = solNum + 1;

                                                solutionClt[solutionClt.size()-1].append(dpndt_var  ==  (temit)*exp( I*tw_coordiPhaseSubs ));
                                            }
                                            else // for handling solutions with conditions
                                            {
                                                solutions<<"solution #"<<solNum<<"  " << dpndt_var << " = "  << (temit)*exp( I*tw_coordiPhaseSubs ) <<" (with condition(s) ";

                                                lst temList = {};
                                                temList.append(dpndt_var  ==  (temit)*exp( I*tw_coordiPhaseSubs ));
                                                for(size_t i=1; i<nops(*it2); i++)
                                                {
                                                    solutions<<(*it2).op(i)<<", ";
                                                    temList.append((*it).op(i));
                                                }

                                                solutions<<");"<<endl;
                                                solNum = solNum + 1;

                                                solutionClt[solutionClt.size()-1].append(temList);

                                            }
                                        }

                                    }
                                }
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
                    //catch(std::domain_error){solutions<<"std::domain_error" <<endl;}


                }




                solutions << endl ;

            }
        }




    }

    solutions <<  "\n" << out << endl;
    solutions << out << "\n" << endl;


/////////////////////////////////////////

    const std::chrono::time_point<std::chrono::system_clock> endTime = chrono::high_resolution_clock::now();
    auto dur = endTime-beginTime;
    cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;
    solutions << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;

    #ifdef GiNaCDE_gui
    stringstream temstr;

    temstr << "The results are written in " << filename <<" file.             " << "Time:"
    << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" ;

    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, &temstr.str()[0]);
    //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,&er);
    /**if(er!=NULL)
    {
        cout << er->message<<endl;
    }**/

    gtk_widget_set_sensitive(window,TRUE);

    #endif // GiNaCDE_gui

    return 0;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////










