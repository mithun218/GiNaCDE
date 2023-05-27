
/** @file utility.cpp
 *
 *  Implementation of some GiNaCDE's usefull utilities used in other files. */



#include <ginac/ginac.h>
#include <math.h>
#include <sstream>
#include<cln/exception.h>
#include "utility.h"
#include "outform.h"
#include "inifcns.h"
#include "simplify.h"
#include "desolve.h"

using namespace std;
using namespace GiNaC;

const symbol symb_=symbol("symb_"),factSymb_ = symbol("factSymb_") ,_FAIL=symbol("_FAIL");

const ex _ex0 = (ex)0, _ex1 = (ex)1, _ex2 = (ex)2, _ex_1 = (ex)-1, _ex_2 = (ex)-2,
 _ex_1_2 = (ex)-1/2, _ex1_2 = (ex)1/2, _ex1_4 = (ex)1/4, _ex_1_4 = (ex)-1/4;

std::chrono::time_point<std::chrono::system_clock> beginTime;

string CurrentPath, filename = "outputResults.txt";


exmapexvec var_depend;

std::vector<lst> solutionClt; lst constraints = {};

symbol_finderc symbol_finder;
dependc depend;

bool has_only_digits(const string s)
{
    return s.find_first_not_of("0123456789") == string::npos;
}
/////////////////////////////////////////////////////////////

vector<string> split (const string &s, char delim)
{
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}
/////////////////////////////////////////////////////////////

/** converting number containing decimal point to rational, 1.252 = 1252/1000 **/
ex DoubleToRational(ex _a)
{
    ostringstream s;
    s<<_a;
    string tem1 = s.str();

    if((tem1.find('E')) != std::string::npos)
    {
        vector<string> ve = split (tem1, 'E');
        size_t zeronum=stoi(ve[1]);
        if(ve[0].find('.') != std::string::npos)
        {
            vector<string> vd = split (ve[0], '.');
            const size_t decipos = vd[1].length();
            zeronum=zeronum-decipos;
            ve[0]=replacestring(ve[0],".","");
        }

        string strCreate;
        if(zeronum>0)
            strCreate = ve[0]+"*10^"+to_string(zeronum);
        else
            strCreate = ve[0]+"*10^(-"+to_string(zeronum)+")";

        return reader(strCreate);
    }

    int j = 0; // position of decimal point
    bool cond = true;
    auto itr = tem1.rbegin();
    do
    {
        if(*itr == '.')
        {
            cond = false;
        }
        else
        {
            j = j + 1;
            itr++;
        }
    }while(cond && itr != tem1.rend());

   if(!cond)
   {
        tem1.erase(find(tem1.begin(),tem1.end(), '.'));
        return reader(tem1)/pow(_ex1*10,j);
   }
   return _a;
}


/** converting number containing decimal point to rational, 1.252 = 1252/1000 **/
ex dorat::operator()(const ex& _e)
{
    if (is_a<numeric>(_e) &&
        !(_e).info(info_flags::rational) && !(_e).info(info_flags::crational))
    {
        israt = false;

        if((_e.info(info_flags::positive) && _e < Gtolerance) ||
                                              (_e.info(info_flags::negative) && -_e < Gtolerance))
                return _ex0;

        if((_e).info(info_flags::real))
            return DoubleToRational(_e);
        else
        {
            return DoubleToRational(real_part(_e))
                   + I*DoubleToRational(imag_part(_e));
        }
    }

    return _e.map(*this);
}


/** replacing power having "add" base with generated symbols **/ // used in Factor
ex powBaseSubs::operator()(const ex& _e)
{
    if((is_a<power>(_e)&&is_a<numeric>(_e.op(1))&&(is_a<add>(_e.op(0))))
            ||(isNu &&is_a<power>(_e)&&_e.op(1).info(info_flags::real) && _e.op(1) < _ex0 && is_a<add>(_e.op(0))))
    {
        if((!exprToSymMap.empty() && exprToSymMap.find(_e.op(0))==exprToSymMap.end())
           ||exprToSymMap.empty())
        {
            j=j+1;

            str = "powSubs_" + to_string(j);
            expr = reader(str);
            exprToSymMap[_e.op(0)]=expr;
        }

        if(!exprToSymMap.empty() && exprToSymMap.find(_e.op(0))!=exprToSymMap.end())
            return pow(exprToSymMap[_e.op(0)],_e.op(1));
    }

    else if(!isNu&&is_a<add>(_e))
        addNum=(addNum+nops(_e))-1;

    return _e.map(*this);
};


/** Polynomial factorization with number containing decimal point and higher degree with "add" base.
    prefactor terms are skipping. **/
ex Factor(const ex& expr)
{
    ex temex;
    if(expr.has(I)) // avoiding I in function. Inbuilt factor does not work when I is present.
    {
        replaceI replaceI;
        temex = replaceI(expr);
    }
    else
        temex = expr;

    if( conjuFreee(temex) == temex && conjuFreee(conjugate(temex)) == temex && !is_a<numeric>(temex) ) // avoiding factor containing conjugate and numerics
    {
        powBaseSubs powSubs(0); //replacing power having "add" base with generated symbols
        fracPowBasSubsFactor fracPow; // replacing base having fractional power.
        fracPow.set();
        temex=fracPow(temex);

        temex = powSubs(temex);// avoiding higher degree with "add" base

        //dorat dorat; // avoiding double form
       // dorat.set();

        try
        {
            if(powSubs.addNum<=addNumFrFactr)
            {
                //cout<<"temexpr_6 "<< temex<<endl;
                temex = factor((temex));
                //temex = collect_common_factors((temex));
                //cout<<"temexpr_7"<< temex<<endl;
            }
            else
                temex = collect_common_factors((temex));
        }
        catch (cln::runtime_exception){} // GiNaC has a bug in factor(). Sometime this error occur.
        catch (std::invalid_argument){}


        temex = genSymbSubs(temex,powSubs.exprToSymMap);
        if(!fracPow.baseClt.empty())
        {
            temex = genSymbSubs(temex,fracPow.baseClt);
        }

        if(temex.has(symb_))
            temex = temex.subs({symb_ == I});

        //if(dorat.israt)
        return temex;
       // else
       //   return evalf(temex);


    }
    else
        return expr;
}

/** calculating gcd of list of expressions (also support rational number: such as {8/7,b/7*b} ) **/
ex Gcd(lst exp)
{
    bool err = false; //isr = true;
    ex tem;
    //dorat dorat; // avoiding double form
    //dorat.set();
    //exp[0] = dorat(exp[0]);
    //if(!dorat.israt)
    //    isr = false;
    while(nops(exp)>=2 && !err)
    {
        try
        {
            //exp[1] = (exp[1]);
            //if(!dorat.israt)
            //    isr = false;
            tem = (gcd(numer(exp[0]),numer(exp[1])))/(lcm(denom(exp[0]),denom(exp[1])));
            exp.remove_first();
            exp.remove_first();
            exp.prepend(tem);
        }
        catch(std::invalid_argument){err = true; exp[0]= _FAIL;}

    }
    //if(isr)
        return exp[0];
    //else
    //    return evalf(exp[0]);
}

/** calculating lcm of list of expressions **/
ex Lcm(lst exp)
{
    bool err = false; //isr = true;
    ex tem;
    //dorat dorat; // avoiding double form
    //dorat.set();
    //exp[0] = dorat(exp[0]);
    //if(!dorat.israt)
    //   isr = false;
    while(nops(exp)>=2 && !err)
    {
        try
        {
            //exp[1] = dorat(exp[1]);
            //if(!dorat.israt)
             //   isr = false;
            tem = (lcm(numer(exp[0]),numer(exp[1])))/(gcd(denom(exp[0]),denom(exp[1])));
            exp.remove_first();
            exp.remove_first();
            exp.prepend(tem);
        }
        catch(std::invalid_argument){err = true; exp[0]= _FAIL;}

    }
    //if(isr)
        return exp[0];
    //else
     //   return evalf(exp[0]);
}


/** replacing I by _symb **/
ex replaceI::operator()(const ex& _e)
{
    if (is_a<numeric>(_e) && !_e.info(info_flags::real))
    {
       return symb_*_e.imag_part() + _e.real_part();
    }
    return _e.map(*this);
}
//////////////////////////////////////////////////////////////////////////////////

int dependc::operator()(const ex& sym, const lst& expr_)
{
   exvector vex1;
   for (auto it = expr_.begin(); it != expr_.end(); it++)
   {
       if( !vex1.empty() )
       {
            if(find(vex1.begin(), vex1.end(), *it) == vex1.end())
                vex1.push_back(*it);
       }
       else
       {
           vex1.push_back(*it);
       }
   }

   for(auto it = var_depend.begin(); it != var_depend.end(); it++ )
   {
      if (it->first == sym)
      {
          exvector vex2 = var_depend[sym], vex3;
          //vex2.insert(vex2.end(), vex1.begin(),vex1.end());
          copy_if( vex1.begin(), vex1.end(), back_inserter(vex3), [&vex2](const ex& arg ){return (find(vex2.begin(), vex2.end(), arg) == vex2.end());});
          if(!vex3.empty())
            vex2.insert(vex2.end(), vex3.begin(),vex3.end());
          var_depend[sym] = vex2;
          return 0;
       }
   }
   var_depend[sym] = vex1;

    return 0;
}
//////////////////////////////////////////////////////////////////////////////////

ex symbol_finderc::operator()(const ex& _e)
{
    if (is_a<symbol>(_e))
    {
       symbols.insert(_e);
       return _e.map(*this);
    }

    return _e.map(*this);
}
//////////////////////////////////////////////////////////////////////////////////

exset symbols(const ex& expr_)
{
    symbol_finder.clear();
    symbol_finder(expr_);
    return symbol_finder.symbols;
}
//////////////////////////////////////////////////////////////////////////////////
/** collecting power of each base from pow argument, excludes similar power**/
ex basepow_clt::operator()(const ex& _e)
{
    if (is_a<power>(_e) &&
        _e.op(0).has(_var))
    {
        lst tem = basepow[_e.op(0)];
        if(nops(tem))
        {
            const exset temset = exprchk;
            exprchk.insert(_e.op(1));
            if(temset.size() != exprchk.size())
                basepow[_e.op(0)] = tem.append(_e.op(1));
        }
        else
            basepow[_e.op(0)] = (lst){_e.op(1)};
        return _e;
    }
    else if(is_a<symbol>(_e) &&
        _e==_var)
    {
        lst tem = basepow[_e];
        if(nops(tem))
        {
            const exset temset = exprchk;
            exprchk.insert(_ex1);
            if(temset.size() != exprchk.size())
                basepow[_e] = tem.append(_ex1);
        }
        else
            basepow[_e] = (lst){_e};
        return _e;
    }


    return _e.map(*this);
}

/** Checking polynomial with variable **/
ex polycheck::operator()(const ex& _e)
{
    if(is_a<power>(_e) && _e.op(0).has(_var) &&
       (((_e.op(1)).info(info_flags::numeric) && !(_e).op(1).info(info_flags::posint))
        || !_e.op(1).info(info_flags::numeric)))
        polytype = false;
    if(is_a<GiNaC::function>(_e) && (_e.op(0).has(_var) || (nops(_e) > 1 && _e.op(1).has(_var))) )
        polytype = false;

    return _e.map(*this);
}

/** Checking presence and type of functions with given variable and checking conjugate in equation **/
ex funcpresent::operator()(const ex& _e)
{
    if(is_a<GiNaC::function>(_e) && nops(_e) == 1)
    {
        if(_e.op(0).has(_var))
        {
            funcwtvar.insert(_e);
            funcpresence = true;
        }
        else if((_e.op(0)) == conjugate(_e)) // avoiding conjugate as it gives error when factor
            funcpresence = true;
        else
            func.insert(_e);
    }
    else if(is_ex_the_function(_e, Diff))
    {
        exmap repls;
        match(_e,Diff(wild(0),wild(1),wild(2)),repls);
        exset symbClt = symbols(repls[wild(0)]);
        for(auto itr = symbClt.begin(); itr != symbClt.end(); itr++)
        {
            exvector dpndClt = depend.get(*itr);
            if(find(dpndClt.begin(), dpndClt.end(), _var) != dpndClt.end())
                funcpresence = true;
        }
     }

     if(is_a<power>(_e))
     {
         if((_e.op(1)).has(_var)) // avoiding solution variable in power
            varInPow = true;
     }


    return _e.map(*this);
}

/// doing conjugate free ///
ex conjuFree::operator()(const ex& _e)
{
    if( nops(_e) == 1 && _e.op(0) == conjugate(_e) )
        return _e.op(0);
    else
        return _e.map(*this);
}

conjuFree conjuFreee;


/// Checking polynomial with variable in fim (avoiding functions other than Diff)///
ex polycheckFim::operator()(const ex& _e)
{
    if(is_a<power>(_e) && _e.op(0).has(_var) && (_e.op(1)).info(info_flags::numeric)
      && !(_e).op(1).info(info_flags::posint))
        polytype = false;
    if(is_a<power>(_e) && _e.op(0).has(_var) && !(_e.op(1)).info(info_flags::numeric))
        polytype = false;
    if(is_a<GiNaC::function>(_e) && _e.op(0).has(_var) && !is_ex_the_function(_e, Diff) && !((_e.op(0)) == conjugate(_e)))
        polytype = false;

    return _e.map(*this);
}

/// Collecting power of a variable in fim (this function has been used for auto-evaluation of NValue in desolve.cpp )///
ex powClt::operator()(const ex& _e)
{
    if(is_a<power>(_e) && _e.op(0) == _var)
        powers.push_back(_e.op(1));
    return _e.map(*this);
}

/** this function substitute generated symbols from exmap **/
ex genSymbSubs(const ex& _e, const exmap& highDegSubsClt)
{
    ex _y = _e;
    if(!highDegSubsClt.empty())
    {
        exmap reversed;
        for(auto itr=highDegSubsClt.begin();itr!=highDegSubsClt.end();itr++)
            reversed[itr->second]=(itr->first);

       _y=(_y.subs(reversed));
    }

    return _y;
}

/** Getting lst of coefficients from all terms where present _var.
isCltPowZero = true allow to get coefficients of _var^0. **/
lst collectAllCoeff(const ex& _expr, const lst& _var, const bool& isCltPowZero, exmap& _varsWtcoeff)
{
    ex coeffClt = _ex1, varsClt = _ex1;
    exmap varsWtcoeff;
    lst totalCoeffClt={};
    if(is_a<add>(_expr))
    {

        for(unsigned i = 0; i < nops(_expr); i++) // add
        {
            if(is_a<mul>(_expr.op(i)))
            {
                for(unsigned j = 0; j < nops(_expr.op(i)); j++)
                {
                    unsigned k = 0;
                    do
                    {
                        if(!_expr.op(i).op(j).has(_var.op(k)))
                        {
                            coeffClt = coeffClt*_expr.op(i).op(j);
                        }
                        else if(_expr.op(i).op(j).has(_var.op(k)))
                        {
                            varsClt = varsClt*_expr.op(i).op(j);
                        }

                        k++;

                    }while(k<nops(_var));
                 }
            }
            else
            {
                unsigned k = 0;
                do
                {
                    if(!_expr.op(i).has(_var.op(k)))
                    {
                        coeffClt = coeffClt*_expr.op(i);
                    }
                    else if(_expr.op(i).has(_var.op(k)))
                    {
                        varsClt = varsClt*_expr.op(i);
                    }

                    k++;

                }while(k<nops(_var));
            }

            if(!isCltPowZero && varsClt != _ex1)
                varsWtcoeff[varsClt] = varsWtcoeff[varsClt]+coeffClt;
            else if(isCltPowZero)
                varsWtcoeff[varsClt] = varsWtcoeff[varsClt]+coeffClt;
            coeffClt=_ex1;
            varsClt=_ex1;
        }


    }
    else
    {
        if(is_a<mul>(_expr))
        {
            for(unsigned j = 0; j < nops(_expr); j++)
            {
                unsigned k = 0;
                do
                {
                    if(!_expr.op(j).has(_var.op(k)))
                    {
                        coeffClt = coeffClt*_expr.op(j);
                    }
                    else if(_expr.op(j).has(_var.op(k)))
                    {
                        varsClt = varsClt*_expr.op(j);
                    }

                    k++;

                }while(k<nops(_var));
             }
        }
        else
        {
            unsigned k = 0;
            do
            {
                if(!_expr.has(_var.op(k)))
                {
                    coeffClt = coeffClt*_expr;
                }
                else if(_expr.has(_var.op(k)))
                {
                    varsClt = varsClt*_expr;
                }

                k++;

            }while(k<nops(_var));
        }

        if(!isCltPowZero && varsClt != _ex1)
            varsWtcoeff[varsClt] = varsWtcoeff[varsClt]+coeffClt;
        else if(isCltPowZero)
            varsWtcoeff[varsClt] = varsWtcoeff[varsClt]+coeffClt;
        coeffClt=_ex1;
        varsClt=_ex1;

    }

    if(!varsWtcoeff.empty())
    {
        _varsWtcoeff=varsWtcoeff;

        for(auto itr = varsWtcoeff.begin(); itr != varsWtcoeff.end(); itr++)
        {
            totalCoeffClt.append(itr->second);
        }
    }

    return totalCoeffClt;
}


/** Getting numerator and denominator. **/
ex Numer_Denom(const ex& _expr)
{
    fracPowBasSubsE.set();
    const ex temexpr_=fracPowBasSubsE(_expr);
    ex nude;
    lst nudeClt;
    nude = temexpr_.numer_denom();

    if(!fracPowBasSubsE.baseClt.empty())
    {
        nudeClt.append( genSymbSubs(nude.op(0),fracPowBasSubsE.baseClt));
        nudeClt.append(genSymbSubs(nude.op(1),fracPowBasSubsE.baseClt));
        fracPowBasSubsE.set();

        return nudeClt;
    }

    return nude;
}

#ifdef GiNaCDE_gui
void resultsinDialog(stringstream& solutions)
{

    GtkTextBuffer *buffer;
    GtkTextIter start, last;

    GtkWidget *dialog=gtk_dialog_new();
    gtk_window_set_decorated(GTK_WINDOW(dialog),TRUE);
    gtk_window_set_title(GTK_WINDOW(dialog),"Solutions with calculating steps->");
    gtk_window_set_modal( GTK_WINDOW(dialog), TRUE );
    gtk_window_set_resizable(GTK_WINDOW(dialog), TRUE);
    gtk_window_set_transient_for( GTK_WINDOW(dialog), GTK_WINDOW(window));
    gtk_widget_set_size_request(GTK_WIDGET(dialog),1024,700);

    GtkWidget *textview3 = gtk_text_view_new();
    GtkWidget *swin1 = gtk_scrolled_window_new (NULL, NULL);
    gtk_container_set_border_width (GTK_CONTAINER (swin1), 5);
    gtk_container_add (GTK_CONTAINER (swin1), textview3);

    gtk_box_pack_start((GtkBox*)(GtkDialog*)(gtk_dialog_get_content_area(GTK_DIALOG(dialog))),swin1,TRUE,TRUE,0);
    //gtk_widget_set_tooltip_text(textview3,"Please enter the system of equations here.");
    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview3));
    gtk_text_buffer_get_start_iter(buffer, &start);
    gtk_text_buffer_get_end_iter(buffer, &last);
    gtk_text_buffer_insert(buffer, &start, &solutions.str()[0],-1);

    gtk_widget_show_all(dialog);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
    return ;
}
#endif //GiNaCDE_gui
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////




