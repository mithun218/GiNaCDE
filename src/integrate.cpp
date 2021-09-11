
/** @file integrate.cpp
 *
 *  Implementation of some basic integration process. */



#include <ginac/ginac.h>
#include "utility.h"
#include "inifcns.h"
#include "simplify.h"
#include "derivative.h"
#include "integrate.h"

using namespace std;
using namespace GiNaC;

//integratec integrate;
/** finding constant terms in integrand **/
/******************************************/
ex integratec::find_const(const ex e_)
{
    ex expr = _ex1;
    if (is_a<mul>(e_))
    {
        for(unsigned i = 0; i != nops(e_); i++ )
        {
            if (!e_.op(i).has(this->var))
                expr = expr*e_.op(i);
        }
    }
    else
    {
        if (!e_.has(this->var))
            expr = e_;
    }
    return expr;
}

/** initial setup of integrand **/
/********************************/
ex integratec::setup(const ex expr_)
{
    if(!is_a<symbol>(this->var))
        return expr_;
    exset sex;
    exmap mexsymb;
    sex = symbols(expr_);
    for (auto it1 = sex.begin(); it1 != sex.end(); it1++)
    {
        auto it2 = var_depend.find(*it1);
        if (it2 != var_depend.end())
        {
            auto it = find(var_depend[*it1].begin(), var_depend[*it1].end(), this->var);
            if (it != var_depend[*it1].end())
            {
                mexsymb[*it1] = Diff(*it1, this->var, _ex0);
            }
        }
    }

    if(!mexsymb.empty())
        return  expr_.subs(mexsymb, subs_options::algebraic);
    else
        return expr_;
}

/** Integration of some functions **/
/*************************************/
ex integratec::func_inte(const ex expr_)
{
    exmap repls;

    /** polynomial **/
    if (expr_ == var)
        return _ex1_2*pow(var, 2);
    if ((1/expr_).is_polynomial(this->var) && (1/expr_).degree(this->var) == 1 )
        return log(1/expr_)/((1/expr_).coeff(this->var, 1));
    expr_.match(pow(var,wild(0)), repls);
    if (!repls.empty() && !repls[wild(0)].has(this->var))
        return pow(var, repls[wild(0)] + _ex1)/(repls[wild(0)] + _ex1);


    /** Exponential **/
    repls.clear();
    expr_.match(exp(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return expr_/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    /** logarithm **/
    expr_.match(log(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return repls[wild(0)]*(log(repls[wild(0)]) + _ex_1)/((repls[wild(0)]).coeff(this->var, 1));

    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 2 )
    {
         const ex a_ =  (repls[wild(0)]).coeff(this->var, 2), b_  = (repls[wild(0)]).coeff(this->var, 1),
         c_ = (repls[wild(0)]).coeff(this->var, 0);
         if ((4*a_*c_ - pow(b_,2)) > _ex0 )
            return this->var*log(repls[wild(0)])-2*this->var+(_ex1_2)*b_*log(repls[wild(0)])/a_+
            4*atan((2*a_*this->var+b_)/sqrt(4*a_*c_-pow(b_,2)))*c_/sqrt(4*a_*c_-pow(b_,2))-
            atan((2*a_*this->var+b_)/sqrt(4*a_*c_-pow(b_,2)))*pow(b_,2)/(sqrt(4*a_*c_-pow(b_,2))*a_);
    }

    repls.clear();
    /** Derivative **/
    expr_.match(Diff(wild(0), wild(1), wild(2)), repls);
    if (!repls.empty())
    {
        if (repls[wild(2)] > _ex0 && repls[wild(1)] == this->var)
            return Diff(repls[wild(0)], repls[wild(1)], repls[wild(2)] - 1);
        else
            return Integrate(expr_, this->var, _ex1);
    }

    repls.clear();
    /** Integration **/
    expr_.match(Integrate(wild(0), this->var, wild(1)), repls);
    if (!repls.empty())
        return Integrate(repls[wild(0)], this->var, repls[wild(1)] + 1);

    repls.clear();
    /** trigonometric functions **/
    expr_.match(sin(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return -cos(repls[wild(0)])/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(cos(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return sin(repls[wild(0)])/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(tan(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return  _ex1_2*log(pow(tan(repls[wild(0)]), 2) + _ex1)/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(csc(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return  log(tan(_ex1_2*repls[wild(0)]))/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(sec(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return  ( log(tan(_ex1_2*(repls[wild(0)])) + _ex1) - log(tan(_ex1_2*(repls[wild(0)])) + _ex_1))/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(cot(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return  (log(tan(_ex1_2*(repls[wild(0)]))) - log(pow(tan(_ex1_2*(repls[wild(0)])), 2) + _ex1))/((repls[wild(0)]).coeff(this->var, 1));

   /** Hyperbolic functions sinh, cosh, tanh **/
    repls.clear();
    expr_.match(sinh(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return cosh(repls[wild(0)])/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(cosh(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return sinh(repls[wild(0)])/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(tanh(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return log(cosh(repls[wild(0)]))/((repls[wild(0)]).coeff(this->var, 1));

    /** some algebraic functions **/

    /** asin, acos, atan, acsc, asec, acot **/
    repls.clear();
    expr_.match(asin(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return repls[wild(0)]*asin(repls[wild(0)]) + sqrt(-pow(repls[wild(0)],2) + _ex1)/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(acos(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return repls[wild(0)]*asin(repls[wild(0)]) + sqrt(-pow(repls[wild(0)],2) + _ex1)/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(atan(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return repls[wild(0)]*atan(repls[wild(0)])+(_ex_1_2)*log(pow(repls[wild(0)],2)+_ex1)/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(acsc(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return acsc(repls[wild(0)])*repls[wild(0)]+log(repls[wild(0)]+repls[wild(0)]*sqrt(_ex1+_ex_1/pow(repls[wild(0)],2)))/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(asec(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return asec(repls[wild(0)])*repls[wild(0)]-log(repls[wild(0)]+repls[wild(0)]*sqrt(_ex1+_ex_1/pow(repls[wild(0)],2)))/((repls[wild(0)]).coeff(this->var, 1));

    repls.clear();
    expr_.match(acot(wild(0)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1 )
        return acot(repls[wild(0)])*repls[wild(0)]+(_ex1_2)*log(pow(repls[wild(0)],2)+_ex1)/((repls[wild(0)]).coeff(this->var, 1));


    return Integrate(expr_, this->var, _ex1);

}

/** integration by parts method **/
/**************************************/

ex integratec::do_inte(const ex expr_)
{
    if (!is_a<mul>(expr_))
        return this->func_inte(expr_);

    exmap repls;

    //x^n*log(x)
    expr_.match(pow(var,wild(0))*log(var), repls);
    if (!repls.empty() && !repls[wild(0)].has(var) && repls[wild(0)] > _ex0)
        return pow(var,(_ex1 + repls[wild(0)]))*(_ex_1 + (_ex1 + repls[wild(0)])*log(var))/pow((1 +
                repls[wild(0)]),2);

    //exp(x)*sin(x)
    repls.clear();
    expr_.match(exp(wild(0))*sin(wild(1)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1
        && (repls[wild(1)]).is_polynomial(this->var) && (repls[wild(1)]).degree(this->var) == 1)
        return exp(repls[wild(0)])*(sin(repls[wild(1)])*(repls[wild(0)]).coeff(this->var, 1)-(repls[wild(1)]).coeff(this->var, 1)*cos(repls[wild(1)]))/(pow((repls[wild(0)]).coeff(this->var, 1),2)+pow((repls[wild(1)]).coeff(this->var, 1),2));

    //exp(x)*cos(x)
    repls.clear();
    expr_.match(exp(wild(0))*cos(wild(1)), repls);
    if (!repls.empty() && (repls[wild(0)]).is_polynomial(this->var) && (repls[wild(0)]).degree(this->var) == 1
        && (repls[wild(1)]).is_polynomial(this->var) && (repls[wild(1)]).degree(this->var) == 1)
        return exp(repls[wild(0)])*((repls[wild(1)]).coeff(this->var, 1)*sin(repls[wild(1)])+cos(repls[wild(1)])*(repls[wild(0)]).coeff(this->var, 1))/(pow((repls[wild(0)]).coeff(this->var, 1),2)+pow((repls[wild(1)]).coeff(this->var, 1),2));

    //Diff*u
    repls.clear();
    expr_.match(Diff(Diff(wild(0),wild(1),0), wild(1),1)*Diff(wild(0),wild(1),0), repls);
    if(!repls.empty())
    {
        if(repls[wild(1)] != var)
            return Integrate(expr_, this->var, _ex1);
        return pow(repls[wild(0)],2)/(2);
    }

    // Diff*u^n
    repls.clear();
    expr_.match(Diff(Diff(wild(0),wild(1),0),wild(1),1)*pow(Diff(wild(0),wild(1),0),wild(2)), repls);
    if(!repls.empty())
    {
        if(repls[wild(1)] != var)
            return Integrate(expr_, this->var, _ex1);
        return pow(repls[wild(0)],repls[wild(2)]+1)/(repls[wild(2)]+1);
    }

    //Diff*Diff
    repls.clear();
    expr_.match(Diff(Diff(wild(0),wild(1),0), wild(1), wild(2))*Diff(Diff(wild(0),wild(1),0), wild(1), wild(4)), repls);
    expr_.match(Diff(Diff(wild(0),wild(1),0), wild(1),wild(2))*Diff(wild(0),wild(1),0), repls);
    if (!repls.empty())
    {
        //avoiding integration of partial derivaives
        if(repls[wild(1)] != var)
            return Integrate(expr_, this->var, _ex1);

        vector<int> dordr;
        if(!(repls[wild(2)]).info(info_flags::integer))
        {
            dorat dorat; // avoiding double form
            dorat.set();
            if(!dorat(repls[wild(2)]).info(info_flags::integer))
                return _FAIL;
            dordr.push_back(ex_to<numeric>(dorat(repls[wild(2)])).to_int());
        }
        else
            dordr.push_back(ex_to<numeric>(repls[wild(2)]).to_int());

        if(!(repls[wild(4)]).info(info_flags::integer))
        {
            dorat dorat; // avoiding double form
            dorat.set();
            if(!dorat(repls[wild(4)]).info(info_flags::integer))
                return _FAIL;
            dordr.push_back(ex_to<numeric>(dorat(repls[wild(4)])).to_int());
        }
        else
            dordr.push_back(ex_to<numeric>(repls[wild(4)]).to_int());

        if((dordr[0] + dordr[1])%2 == 0)
            return Integrate(expr_, this->var, _ex1);

        sort(dordr.begin(), dordr.end());

        exmap rep, rep1;
        ex exmax, exmin;
        expr_.op(0).match(Diff(wild(5), wild(6), wild(7)), rep);

        int wild7int;
        if(!(rep[wild(7)]).info(info_flags::integer))
        {
            dorat dorat; // avoiding double form
            dorat.set();
            if(!dorat(rep[wild(7)]).info(info_flags::integer))
                return _FAIL;
            wild7int=ex_to<numeric>(dorat(rep[wild(7)])).to_int();
        }
        else
            wild7int=ex_to<numeric>(rep[wild(7)]).to_int();

        if(wild7int == dordr[0])
        {
            exmin = expr_.op(0);
            exmax = expr_.op(1);
        }
        else
        {
            exmin = expr_.op(1);
            exmax = expr_.op(0);
        }

        rep.clear();

        exmin.match(Diff(wild(5), wild(6), wild(7)), rep);
        exmax.match(Diff(wild(5), wild(6), wild(7)), rep1);

        ex expr = _ex0;
        int imax = (dordr[1] - 1 - dordr[0])/2;

        for(int i = 0; i <= imax; i++)
        {
            if(i == imax)
            {
                expr = expr + _ex1_2*pow(_ex_1, i)*Diff(rep[wild(5)], rep[wild(6)], rep[wild(7)] + i)*
                    Diff(rep1[wild(5)], rep1[wild(6)], rep1[wild(7)] - 1 - i);
            }
            else
            {
                expr = expr + pow(_ex_1, i)*Diff(rep[wild(5)], rep[wild(6)], rep[wild(7)] + i)*
                    Diff(rep1[wild(5)], rep1[wild(6)], rep1[wild(7)] - 1 - i);
            }

        }
        rep.clear();
        rep1.clear();

        return expr;
    }


    repls.clear();
    //with x^n*Diff(wild(0)) funtions
    expr_.match(wild(0)*Diff(wild(1), wild(2), wild(3)), repls);
    //x^n*sin(a*x+b)
    expr_.match(wild(0)*sin(wild(1)), repls);
    //x^n*cos(a*x+b)
    expr_.match(wild(0)*cos(wild(1)), repls);
    //x^n*exp(a*x+b)
    expr_.match(wild(0)*exp(wild(1)), repls);

    if (!repls.empty())
    {
        //avoiding integration of partial derivaives
        if(expr_.has(Diff(wild(1), wild(2), wild(3))) && repls[wild(2)] != var)
            return Integrate(expr_, this->var, _ex1);

        ex funcexpr;

        funcexpr = Simplify(expand(expr_/repls[wild(0)]));

        if( !funcexpr.has(Diff(wild(0), wild(1), wild(2))) )
        {
            if (!repls[wild(0)].is_polynomial(this->var) || !repls[wild(1)].is_polynomial(this->var)
                || !((repls[wild(1)]).degree(this->var) == 1))
                return Integrate(expr_, this->var, _ex1);
        }
        else
        {
            if (!repls[wild(0)].is_polynomial(this->var))
                return Integrate(expr_, this->var, _ex1);
        }

        ex expr1, expr2, expr = _ex0, constterm;
        exvector exvec1, exvec2;
        exvec1.push_back(repls[wild(0)]);
        expr1 = repls[wild(0)];
        expr2 = funcexpr;
        do
        {
            constterm = this->find_const(expr2);
            expr2 = constterm*this->func_inte( Simplify(expand(expr2/constterm )));
            if (expr2.has(Integrate(wild(0), wild(1), wild(2))) && exvec1.size() == 1)
                return Integrate(expr_, this->var, _ex1 );
            expr1 = pdiff(expr1, this->var, _ex1);

            exvec1.push_back(expr1);
            exvec2.push_back(expr2);
        }while(expr1 != _ex0 && !expr2.has(Integrate(wild(0), wild(1), wild(2))));

        auto it1 = exvec1.begin();
        auto it2 = exvec2.begin();

        for (size_t i = 0; i < exvec1.size() - 1; i++)
        {
            if (!((*next(it2, i)).has(Integrate(wild(0), wild(1), wild(2)))))
                expr = expr + pow(_ex_1, i)*(*next(it1, i))*(*next(it2, i));
            else
                expr = expr + pow(_ex_1, i)*Integrate((*next(it1, i))*(*next(it2, i)), this->var, _ex1);
        }
        return expr;
    }



    return Integrate(expr_, this->var, _ex1);
    //with Integrate
}

/** implementation of the class integratec **/
/**************************************/

ex integratec::operator()(const ex & expr_)
{
    if(expr_ == _ex0)
        return expr_;

    ex expr, constterms;

    expr = expr_;
    expr = Simplify(expand(expr), TrigCombine);
    expr = Simplify(expand(expr), logSimp);
    //expr = evaluate(expr);

    expr = this->setup(expr);
    if (is_a<add>(expr))
    {
        exvector ev;
        for(size_t i=0;i<expr.nops();++i)
        {
            constterms = this->find_const(expr.op(i));
            if (expr.op(i).is_equal(constterms))
                 ev.push_back(expr.op(i)*this->var);
            else
            {
                ev.push_back(Simplify(expand(this->do_inte(Simplify(expand(expr.op(i)/constterms)))*constterms)));
            }
        }
        ex tem = add(ev);
        return zero_order_rem(tem);
    }
    else
    {
        constterms = this->find_const(expr);
        if (expr.is_equal(constterms))
        {
            expr = expr*(this->var);
            return zero_order_rem(expr);
        }
        expr = Simplify(expand(expr/constterms));
        expr = this->do_inte(expr);
        expr = expr*constterms;
        expr = Simplify(expand(expr));
        expr = zero_order_rem(expr);
        return expr;
    }
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////























