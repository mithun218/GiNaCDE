
/** @file solve.cpp
 *
 *  Solving a system of linear or nonlinear polynomial equations (implementation). */



#include <iostream>
#include <bits/stdc++.h>
#include <ginac/ginac.h>
#include<cln/exception.h>
#include "simplify.h"
#include "utility.h"
#include "solve.h"
#include "inifcns.h"

using namespace std;
using namespace GiNaC;


/** counting total number of add in an equation**/
ex totalAddInEq::operator()(const ex &_e)
{
    if(is_a<add>(_e))
        addNum=addNum+(nops(_e)-1);
    return _e.map(*this);
}

/** replacing the "pow" terms where present _var with created symbols **/
ex powBaseSubsWtVar::operator()(const ex& _e)
{
    if(exprToSymMapWtVar.size()>1)
        return _e;

    else if(is_a<power>(_e)&&(_e.op(0).has(_var))) // substituting power which have _var
    {
        if((!exprToSymMapWtVar.empty() && exprToSymMapWtVar.find(_e)==exprToSymMapWtVar.end())
            ||exprToSymMapWtVar.empty())
        {
            j=j+1;

            str = "powBaseSubsWtVar_" + to_string(j);
            expr = reader(str);
            exprToSymMapWtVar[_e]=expr;
        }
        if(!exprToSymMapWtVar.empty() && exprToSymMapWtVar.find(_e)!=exprToSymMapWtVar.end())
            return exprToSymMapWtVar[_e];
    }

    return _e.map(*this);
};

/** replacing the "pow" terms without _var with created symbols. **/
ex powBaseSubsWtoutVar::operator()(const ex& _e)
{
    if(is_a<power>(_e)&&!(_e.op(0).has(_var))) // substituting power which does not have _var
    {
        if((!exprToSymMapWtoutVar.empty() && exprToSymMapWtoutVar.find(_e)==exprToSymMapWtoutVar.end())
            ||exprToSymMapWtoutVar.empty())
        {
            j=j+1;

            str = "powBaseSubsWtout_" + to_string(j);
            expr = reader(str);
            exprToSymMapWtoutVar[_e]=expr;
        }
        if(!exprToSymMapWtoutVar.empty() && exprToSymMapWtoutVar.find(_e)!=exprToSymMapWtoutVar.end())
            return exprToSymMapWtoutVar[_e];
    }

    return _e.map(*this);
};

/** Square root of complex numbers. **/
ex solvec::sqroot(const ex& _exp)
{
    if(!is_a<numeric>(_exp))
        return _ex_1;
    if(_exp.info(info_flags::real))
        return sqrt(_exp);
    else
    {
        ex re, im;
        re = real_part(evalf(_exp));
        im = imag_part(evalf(_exp));
        ex tem = evalf(sqrt(_ex1_2*(sqrt(re*re+im*im) - re)));
        return evalf(_ex1_2*im/tem) + sqrt(_ex_1)*tem;
    }
}


/** Solutions of linear and quadratic polynomial equation **/
lst solvec::quadratic(const ex & equ_, const ex& var_)
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
        ex tem = Simplify(_b*_b - 4*_a*_c);
        if(tem == _ex0)
        {
            quasolu.append(var_ == simplify((_ex1_2*(-_b))/_a));
        }
        else
        {
            quasolu.append(var_ == simplify(Simplify(_ex1_2*(-_b + sqrt(tem)))/_a));
            quasolu.append(var_ == simplify(Simplify(_ex1_2*(-_b - sqrt(tem)))/_a));
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

/** algebraic solutions of cubic polynomial equation by cardan's method **/
lst solvec::cubics(const ex & equ_, const ex & var_)
{
    lst cubicsolu;
    try
    {

        const ex _a = Simplify(equ_.coeff(var_, 3)),
                 _b = Simplify((equ_.coeff(var_, 2))),
                 _c = Simplify((equ_.coeff(var_, 1))),
                 _d = Simplify((equ_.coeff(var_, 0)));

        const ex S=pow((9*_a*_b*_c-27*_a*_a*_d-2*_b*_b*_b)/(54*_a*_a*_a) + sqrt(pow((3*_a*_c-_b*_b)/(9*_a*_a),3) + pow((9*_a*_b*_c-27*_a*_a*_d-2*_b*_b*_b)/(54*_a*_a*_a),2)) ,_ex1/3),
                 T=pow((9*_a*_b*_c-27*_a*_a*_d-2*_b*_b*_b)/(54*_a*_a*_a) - sqrt(pow((3*_a*_c-_b*_b)/(9*_a*_a),3) + pow((9*_a*_b*_c-27*_a*_a*_d-2*_b*_b*_b)/(54*_a*_a*_a),2)) ,_ex1/3);

        cubicsolu.append(var_== Simplify((S+T)-_b/(3*_a)));
        cubicsolu.append(var_== Simplify(-(S+T)/2-_b/(3*_a)+I*(sqrt(_ex1*3)/2)*(S-T)));
        cubicsolu.append(var_== Simplify(-(S+T)/2-_b/(3*_a)-I*(sqrt(_ex1*3)/2)*(S-T)));

        for(unsigned i = 0; i < nops(cubicsolu); i++)
        {
            if(cubicsolu.op(i).rhs().info(info_flags::positive) && cubicsolu.op(i).rhs() < Gtolerance)
                cubicsolu[i] = {var_ == _ex0};
            else if(cubicsolu.op(i).rhs().info(info_flags::negative) && -cubicsolu.op(i).rhs() < Gtolerance)
                cubicsolu[i] = {var_ == _ex0};
        }
    }
    catch(cln::runtime_exception){return cubicsolu;}


    return cubicsolu;
}
/**Numerical solutions of cubic polynomial equation **/
lst solvec::Ncubics(const ex & equ_, const ex& var_)
{
    lst cubicsolu;

    const ex _p = equ_.coeff(var_, 3),
            p = Simplify(equ_.coeff(var_, 2)/_p),
            q = Simplify(equ_.coeff(var_, 1)/_p),
            r = Simplify(equ_.coeff(var_, 0)/_p);

    if(!_p.info(info_flags::numeric) || !p.info(info_flags::numeric) || !q.info(info_flags::numeric)
       || !r.info(info_flags::numeric) )
      {
          return {};
      }

    const ex _ex1_27 = (ex)1/27, _ex1_3 = (ex)1/3;

    const ex a = Simplify(_ex1_3*(3*q - p*p)),
             b = Simplify(_ex1_27*(2*p*p*p - 9*p*q + 27*r));

    ex discriminant = Simplify(_ex1_4*b*b + _ex1_27*a*a*a);
    if(is_a<numeric>(b) && is_a<numeric>(discriminant))
    {
        if(discriminant == _ex0)
        {
            if(b > _ex0)
            {
                cubicsolu.append(var_ == -2*sqrt(-_ex1_3*a) - _ex1_3*p );
                cubicsolu.append(var_ == sqrt(-_ex1_3*a) - _ex1_3*p );
                cubicsolu.append(var_ == sqrt(-_ex1_3*a) - _ex1_3*p );
            }
            else if(b < _ex0)
            {
                cubicsolu.append(var_ == 2*sqrt(-_ex1_3*a) - _ex1_3*p );
                cubicsolu.append(var_ == -sqrt(-_ex1_3*a) - _ex1_3*p );
                cubicsolu.append(var_ == -sqrt(-_ex1_3*a) - _ex1_3*p );
            }
            else
            {
                cubicsolu.append(var_ == _ex0 - _ex1_3*p );
                cubicsolu.append(var_ == _ex0 - _ex1_3*p );
                cubicsolu.append(var_ == _ex0 - _ex1_3*p );
            }
        }

        else if(discriminant < _ex0)
        {
            ex phi;

            if(b >= _ex0)
            {
               phi = acos(-sqrt(_ex_1_4*(27*b*b)/(a*a*a)));
               if(is_ex_the_function(phi, acos))
                    phi = acos(evalf(-sqrt(_ex_1_4*(27*b*b)/(a*a*a))));
            }
            else
            {
                phi = acos(sqrt(_ex_1_4*(27*b*b)/(a*a*a)));
               if(is_ex_the_function(phi, acos))
                    phi = acos(evalf(sqrt(_ex_1_4*(27*b*b)/(a*a*a))));
            }

            ex tem = cos(_ex1_3*phi);
            if(is_ex_the_function(tem, cos))
                tem = cos(evalf(_ex1_3*phi));

            cubicsolu.append(var_ == evalf(Simplify(-_ex1_3*p + 2*sqrt(-_ex1_3*a)*(tem))));

            tem = cos(_ex1_3*phi + 2*_ex1_3*(Pi));
            if(is_ex_the_function(tem, cos))
                tem = cos(evalf(_ex1_3*phi + 2*_ex1_3*(Pi)));
            cubicsolu.append(var_ == evalf(Simplify(-_ex1_3*p + 2*sqrt(-_ex1_3*a)*(tem))));

            tem = cos(_ex1_3*phi + 4*_ex1_3*(Pi));
            if(is_ex_the_function(tem, cos))
                tem = cos(evalf(_ex1_3*phi + 4*_ex1_3*(Pi)));
            cubicsolu.append(var_ == evalf(Simplify(-_ex1_3*p + 2*sqrt(-_ex1_3*a)*(tem))));

        }

        else
        {
            ex A, B, tem;

            tem = evalf(Simplify(_ex_1_2*b + sqrt(discriminant)));
            if(tem.info(info_flags::negative))
                A = -evalf(pow((-tem), _ex1_3));
            else
                A = evalf(pow((tem), _ex1_3));
            tem = evalf(Simplify(_ex_1_2*b - sqrt(discriminant)));

            if(tem.info(info_flags::negative))
                B = -evalf(pow((-tem), _ex1_3));
            else
                B = evalf(pow((tem), _ex1_3));

            ex three = 3;
            cubicsolu.append(var_ == evalf(Simplify(A + B - _ex1_3*p)));
            cubicsolu.append(var_ == evalf(Simplify(_ex_1_2*(A + B) + sqrt(_ex_1)*_ex1_2*sqrt(three)*(A - B) - _ex1_3*p)));
            cubicsolu.append(var_ == evalf(Simplify(_ex_1_2*(A + B) - sqrt(_ex_1)*_ex1_2*sqrt(three)*(A - B) - _ex1_3*p)));
        }
    }

    for(unsigned i = 0; i < nops(cubicsolu); i++)
    {
        if(cubicsolu.op(i).rhs().info(info_flags::positive) && cubicsolu.op(i).rhs() < GiNaC::pow(10, -10))
            cubicsolu[i] = {var_ == _ex0};
        else if(cubicsolu.op(i).rhs().info(info_flags::negative) && -cubicsolu.op(i).rhs() < GiNaC::pow(10, -10))
            cubicsolu[i] = {var_ == _ex0};
    }

    return cubicsolu;

}

/**Numerical solutions of quartic polynomial equation (Euler's solution) **/
lst solvec::Nquartic(const ex & equ_, const ex& var_)
{
    lst quarticsolu={};

    const ex a = equ_.coeff(var_, 4),
            b = equ_.coeff(var_, 3)/a,
            c = equ_.coeff(var_, 2)/a,
            d = equ_.coeff(var_, 1)/a,
            e = equ_.coeff(var_, 0)/a;

    if(!a.info(info_flags::numeric) || !b.info(info_flags::numeric) || !c.info(info_flags::numeric)
       || !d.info(info_flags::numeric) || !e.info(info_flags::numeric) )
        {
            return {};
        }

    const ex _ex1_8 = (ex)1/8, _ex1_256 = (ex)1/256, _ex1_16 = (ex)1/16, _ex1_64 = (ex)1/64;

    const ex f = Simplify(c - (3*b*b*_ex1_8)),
            g = Simplify(d + (_ex1_8*b*b*b) - (_ex1_2*b*c)),
            h = Simplify((e) - 3*_ex1_256*pow(b,4) +(_ex1_16*b*b*c) - (_ex1_4*b*d));

    const lst cubicsolu = this->Ncubics(pow(var_, 3) + (_ex1_2*f)*var_*var_ + Simplify(_ex1_16*(f*f - 4*h))*var_ -_ex1_64*g*g, var_);

    if(nops(cubicsolu) == _ex0)
        return quarticsolu;

    lst sqrtli;

    for(auto it =cubicsolu.begin(); it != cubicsolu.end(); it++)
    {
        if(!((*it).info(info_flags::real)))
            sqrtli.append(this->sqroot(*it));
    }
    if(nops(sqrtli) < 2)
    {
       for(auto it = cubicsolu.begin(); it != cubicsolu.end(); it++)
       {
        if((*it) != _ex0)
            sqrtli.append(this->sqroot(*it));
        }
    }

    if(nops(sqrtli) < 2)
     return quarticsolu;

    const ex p = sqrtli.op(0),
            q = sqrtli.op(1),
            r = Simplify(-_ex1_8*g/(p*q)),
            s = Simplify(_ex1_4*b);

    quarticsolu.append(var_ == evalf(Simplify(p + q + r -s)));
    quarticsolu.append(var_ == evalf(Simplify(p - q - r -s)));
    quarticsolu.append(var_ == evalf(Simplify(-p + q - r - s)));
    quarticsolu.append(var_ == evalf(Simplify(-p - q + r - s)));

    for(unsigned i = 0; i < nops(quarticsolu); i++)
    {
        if(quarticsolu.op(i).rhs().info(info_flags::positive) && quarticsolu.op(i).rhs() < Gtolerance)
            quarticsolu[i] = {var_ == _ex0};
        else if(quarticsolu.op(i).rhs().info(info_flags::negative) && -quarticsolu.op(i).rhs() < Gtolerance)
            quarticsolu[i] = {var_ == _ex0};
    }

    return quarticsolu;

}

/** solution of polynomial equation with automatic degree selection **/
lst solvec::polySoluWtAutoDegSelct(const ex& _equ, const ex& _var)
{
    ex temSolu; lst tsolulst;
    ex temex = _equ,temex2;

    bool isContinue = true;
    powBaseSubsWtVar powBaseSubsWtVar(0,_var);

    do
    {
        powBaseSubsWtVar.exprToSymMapWtVar.clear();
        temex = powBaseSubsWtVar(temex);

        if(powBaseSubsWtVar.exprToSymMapWtVar.empty()&&temex.has(_var))
            isContinue = false;
        else if((!powBaseSubsWtVar.exprToSymMapWtVar.empty()&&temex.has(_var))
                ||powBaseSubsWtVar.exprToSymMapWtVar.size()>1)
        {
            temex = genSymbSubs(temex,powBaseSubsWtVar.exprToSymMapWtVar);
            isContinue = false;
        }
        else if(!temex.has(_var)
        && powBaseSubsWtVar.exprToSymMapWtVar.size()==1
        && is_polynomial(temex,powBaseSubsWtVar.exprToSymMapWtVar.begin()->second)
        && degree(temex,powBaseSubsWtVar.exprToSymMapWtVar.begin()->second)==1)
        {
            temSolu = -Simplify(collect_common_factors(coeff(temex,powBaseSubsWtVar.exprToSymMapWtVar.begin()->second,0))/collect_common_factors(coeff(temex,powBaseSubsWtVar.exprToSymMapWtVar.begin()->second,1)));

            if(is_a<power>(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first)
               &is_polynomial(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0),_var)
               &&(degree(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0),_var)==1
               ||degree(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0),_var)==2))
            {

                if(degree(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0),_var)==1)  // Ex: (a*y+b)^2+c, (a*y+b)^3+c
                {
                   lst ttsolu,tttsolu;
                   if(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1).info(info_flags::real)
                    && powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1).info(info_flags::even) )// Even powers return +,- solutions. Ex: (a*y+b)^2+c
                   {
                        ttsolu=quadratic(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0)+pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);
                        if(nops(ttsolu)!=0)
                            tttsolu.append(ttsolu.op(0));
                        ttsolu=quadratic(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0)-pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);
                        if(nops(ttsolu)!=0)
                            tttsolu.append(ttsolu.op(0));

                        return tttsolu;

                   }
                   else // ode powers return + solutions. Ex: (a*y+b)^3+c
                   {
                        ttsolu=quadratic(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0)-pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);
                        if(nops(ttsolu)!=0)
                            tttsolu.append(ttsolu.op(0));

                        return tttsolu;
                   }
                }
                if(degree(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0),_var)==2) // Ex: (a*y^2+b*y++b)^2+c, (a*y^2+b*y++b)^3+c
                {
                   powBaseSubsWtoutVar powBaseSubsWtoutVar(0,_var);
                   temex2 = powBaseSubsWtoutVar(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0));
                   temex2 = expand(temex2);
                   temex2 = genSymbSubs(temex2,powBaseSubsWtoutVar.exprToSymMapWtoutVar);

                   lst ttsolu,tttsolu;
                   if(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1).info(info_flags::real)
                    && powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1).info(info_flags::even) ) // Even powers return +,- solutions. Ex: (a*y^2+b*y++b)^2+c
                   {
                       ttsolu=this->quadratic(temex2+pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);
                       for(unsigned i=0;i<nops(ttsolu);i++)
                            tttsolu.append(ttsolu.op(i));
                       ttsolu=this->quadratic(temex2-pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);
                       for(unsigned i=0;i<nops(ttsolu);i++)
                            tttsolu.append(ttsolu.op(i));

                        return tttsolu;

                   }
                   else// Odd powers return + solutions. Ex: (a*y^2+b*y++b)^3+c
                   {
                   ttsolu=this->quadratic(temex2-pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1))),_var);

                   for(unsigned i=0;i<nops(ttsolu);i++)
                        tttsolu.append(ttsolu.op(i));

                    return tttsolu;
                   }
                }

            }
            else
            {
                temex = powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(0)-pow(temSolu,_ex1/(powBaseSubsWtVar.exprToSymMapWtVar.begin()->first.op(1)));
            }
        }

    }while(isContinue);

    temex = _equ;

    if(is_polynomial(temex,_var)
    && (unsigned)degree(temex,_var)<=expandLevel)
    {
        powBaseSubsWtoutVar powBaseSubsWtoutVar(0,_var);
        temex = powBaseSubsWtoutVar(temex);
        temex = Simplify(expand(temex));

        if(temex==_ex0)
        {
            tsolulst.append(_var==_var);
            return tsolulst;
        }
        else if(!temex.has(_var))
        {
            return tsolulst;
        }

        temex = genSymbSubs(temex,powBaseSubsWtoutVar.exprToSymMapWtoutVar);

        if((temex).degree(_var) < 3 && temex.degree(_var) > 0) // linear and quadratic equation
        {
            tsolulst = this->quadratic((temex), _var);

            if(nops(tsolulst)!=0)
                return tsolulst;
        }
        else if((temex).degree(_var) == 3) // cubic equation
        {
            tsolulst = this->Ncubics((temex), _var);

            if(nops(tsolulst)!=0)
                return tsolulst;
        }
        else if((temex).degree(_var) == 4) // quartic equation
        {
            tsolulst = this->Nquartic((temex), _var);

            if(nops(tsolulst)!=0)
                return tsolulst;
        }
    }

    temex = _equ;

    basepow_clt basepow_clt(_var);
    basepow_clt.clear();
    basepow_clt(temex);

    if(basepow_clt.basepow.size() == 1)  // attempt 2 for getting solutions
    {
        lst temsolu;

        // calculating gcd of collected powers
        lst powlst = basepow_clt.basepow.begin()->second;
        ex GCd = Gcd(powlst);

        if(GCd != _FAIL)
        {

            exmap subsex;
            ex Z_=reader("Z_");
            //for(auto itr = powlst.begin(); itr != powlst.end(); itr++)
            //{
            subsex[(basepow_clt.basepow.begin()->first)] = pow(Z_,Simplify(_ex1/GCd));

            //}
            ex temlow_var_eq;
            try
            {
                temlow_var_eq = Simplify(temex.subs(subsex,subs_options::algebraic));
            }
            catch(GiNaC::pole_error){return tsolulst;}

            if(is_polynomial(temlow_var_eq, Z_) && degree(temlow_var_eq, Z_)<3)
            {
                temsolu=this->quadratic(temlow_var_eq, Z_);

                if(nops(temsolu)!=0)
                {
                    lst ttsolu,tttsolu;
                    for(unsigned i = 0;i<nops(temsolu);i++)
                    {
                        if(is_polynomial(basepow_clt.basepow.begin()->first,_var)
                           &&(degree(basepow_clt.basepow.begin()->first,_var)==1
                           ||degree(basepow_clt.basepow.begin()->first,_var)==2))
                        {

                            if(degree(basepow_clt.basepow.begin()->first,_var)==1)  // Ex: (a*y+b)^2+c, (a*y+b)^3+c
                            {
                               if(GCd.info(info_flags::real)
                                && GCd.info(info_flags::even) )// Ex: (a*y+b)^2+c
                               {
                                    ttsolu=quadratic(basepow_clt.basepow.begin()->first+pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);
                                    if(nops(ttsolu)!=0)
                                        tttsolu.append(ttsolu.op(0));
                                    ttsolu=quadratic(basepow_clt.basepow.begin()->first-pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);
                                    if(nops(ttsolu)!=0)
                                        tttsolu.append(ttsolu.op(0));

                               }
                               else // Ex: (a*y+b)^3+c
                               {
                                    ttsolu=quadratic(basepow_clt.basepow.begin()->first-pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);
                                    if(nops(ttsolu)!=0)
                                        tttsolu.append(ttsolu.op(0));
                               }
                            }
                            if(degree(basepow_clt.basepow.begin()->first,_var)==2) // Ex: (a*y^2+b*y++b)^2+c, (a*y^2+b*y++b)^3+c
                            {
                               powBaseSubsWtoutVar powBaseSubsWtoutVar(0,_var);
                               temex2 = powBaseSubsWtoutVar(basepow_clt.basepow.begin()->first);
                               temex2 = expand(temex2);
                               temex2 = genSymbSubs(temex2,powBaseSubsWtoutVar.exprToSymMapWtoutVar);

                               if(GCd.info(info_flags::real)
                                && GCd.info(info_flags::even) ) // Ex: (a*y^2+b*y++b)^2+c
                               {
                                   ttsolu=this->quadratic(temex2+pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);
                                   for(unsigned i=0;i<nops(ttsolu);i++)
                                        tttsolu.append(ttsolu.op(i));
                                   ttsolu=this->quadratic(temex2-pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);
                                   for(unsigned i=0;i<nops(ttsolu);i++)
                                        tttsolu.append(ttsolu.op(i));

                               }
                               else// Ex: (a*y^2+b*y++b)^3+c
                               {
                               ttsolu=this->quadratic(temex2-pow(temsolu.op(i).rhs(),_ex1/(GCd)),_var);

                               for(unsigned i=0;i<nops(ttsolu);i++)
                                    tttsolu.append(ttsolu.op(i));
                               }
                            }

                        }

                    }

                    if(nops(tttsolu)!=0)
                        return tttsolu;

                }
            }
        }

    }


    return tsolulst;
}


/** getting solutions of a polynomial equation **/
int solvec::one_eq_solutions(const ex& _equ, const ex& _var)
{
    one_eq_solu.remove_all();

    if(!_equ.has(_var))
        return 0;

    lst temSolu;

    if((is_a<power>(_equ)&&is_a<add>(_equ.op(0)))
       ||(is_a<add>(_equ)))
    {
        temSolu = polySoluWtAutoDegSelct(_equ,_var);

        if(nops(temSolu != 0))
        {
            for(auto it = temSolu.begin(); it != temSolu.end(); it++)
            {
                if(!one_eq_solu.has(*it))
                    one_eq_solu.append(*it);
            }
         }
    }
    else
    {
        try
        {

            if((_equ).subs(_var == _ex0)==_ex0)
                one_eq_solu.append(_var == _ex0);
        }
        catch(GiNaC::pole_error){return 0;}
    }

    if(nops(one_eq_solu) == _ex0)
            return 0;
    return 1;
}

/** solution variables are sorted with ascending order of degree of variables in poly equation (eqDivider[j]) **/
lst solvec::varOrderByDeg(const lst& low_var, map<unsigned, ex>& eqDivider, unsigned& eqDividerSz)
{
    lst low_var_sep, low_var_sep_wtotpoly;
    if(eqDividerSz > eqDivider.size())
       return low_var_sep;


    std::vector<std::pair<ex, int>> eq_degnum_pairs;

    do
    {
          for( auto j =  low_var.begin(); j !=  low_var.end(); j++)
          {
             if((eqDivider[eqDividerSz]).has(*j) && is_polynomial(eqDivider[eqDividerSz],*j))
             {
                eq_degnum_pairs.push_back(make_pair(*j, degree(eqDivider[eqDividerSz],*j)));
             }
             else if((eqDivider[eqDividerSz]).has(*j))
                low_var_sep_wtotpoly.append(*j);  // non-polynomial variables are collected
          }
          if(!eq_degnum_pairs.empty())
          {
              sort(eq_degnum_pairs.begin(), eq_degnum_pairs.end(), [=](std::pair<ex, int>& a, std::pair<ex, int>& b){
                  return a.second < b.second;}); // sorting of polynomial variables

              for(auto itr = eq_degnum_pairs.begin(); itr != eq_degnum_pairs.end(); itr++)
              {
                  low_var_sep.append(itr->first);
              }

          }
          if(nops(low_var_sep_wtotpoly) != 0)
          {
              for(auto itr = low_var_sep_wtotpoly.begin(); itr != low_var_sep_wtotpoly.end(); itr++)
              {
                  low_var_sep.append(*itr);
              }
          }
          if(nops(low_var_sep) != 0)
             return low_var_sep;

              eqDividerSz++; // shifting to next equation

    }while(eqDividerSz < eqDivider.size());

   return low_var_sep;
}

// checking present of variables (_var) in _expr
bool solvec::isVarPrsnt(const ex &_expr, const exset &_var)
{

    bool funcPrsnt = false;
    for( auto itr =  _var.begin(); itr !=  _var.end();)
    {

        if(_expr.has(*itr))
        {
            funcPrsnt = true;
            itr = _var.end();
        }
        else
        {
            itr++;
        }
    }

    return funcPrsnt;

}

ex solvec::sysequ_red(const exset& sysequ_, const exset& var_)
{

    /// collecting solutions ///
    if(sysequ_.empty())
    {
        solu.insert(soluClt);
    }
    else
    {
        try
        {
            /// Eliminating variables ///
            // selecting equations with the lowest number of solution variables
           lst low_var, var_in_eq;
           ex low_var_eq, ret, denoma;
           exset sysequ = sysequ_;
           ex subs_simp;
           bool isSolu;

           var_in_eq.remove_all();

           std::vector<size_t> addNumFrEchtrmIneq;
           std::vector<std::pair<ex, unsigned>> maxAddnum;
           std::map<ex,unsigned,ex_is_less> totalAddnum;

           for (auto itr = sysequ.begin(); itr != sysequ.end(); ++itr)
           {
               addNumFrEchtrmIneq.clear();
               bool isVarprsntInEq = false;
               for( auto j =  var_.begin(); j !=  var_.end(); j++) // checking present of minimum one solution-var
               {                                                   // in each equations.
                   if((*itr).has(*j))
                       isVarprsntInEq = true;
               }
               if(!isVarprsntInEq) // if minimum one solution-var not present
                   return _ex_1;   // solutions do not exist.


               if(is_a<add>(*itr))
               {
                   if(isVarPrsnt(*itr, var_))
                   {
                       totalAddInEqV.addNum = 0;
                       totalAddInEqV(*itr);
                       addNumFrEchtrmIneq.push_back(totalAddInEqV.addNum);
                   }
                   else
                   {
                       addNumFrEchtrmIneq.push_back(0);
                   }
               }

               else if(is_a<power>(*itr) && is_a<add>((*itr).op(0)))
               {
                   if(isVarPrsnt((*itr).op(0), var_))
                   {
                       totalAddInEqV.addNum = 0;
                       totalAddInEqV((*itr).op(0));
                       addNumFrEchtrmIneq.push_back(totalAddInEqV.addNum);
                   }
               }

               else if(is_a<mul>(*itr))
               {
                   for (size_t i = 0; i<nops(*itr);i++)
                   {
                       if(is_a<add>((*itr).op(i)))
                       {
                           if(isVarPrsnt((*itr).op(i), var_))
                           {
                               totalAddInEqV.addNum = 0;
                               totalAddInEqV((*itr).op(i));
                               addNumFrEchtrmIneq.push_back(totalAddInEqV.addNum);
                           }
                       }
                       else if(is_a<power>((*itr).op(i)) && is_a<add>((*itr).op(i).op(0)))
                       {
                           if(isVarPrsnt((*itr).op(0), var_))
                           {
                               totalAddInEqV.addNum = 0;
                               totalAddInEqV((*itr).op(i).op(0));
                               addNumFrEchtrmIneq.push_back(totalAddInEqV.addNum);
                           }
                       }
                       else
                       {
                           if(isVarPrsnt((*itr).op(i), var_))
                           {
                               addNumFrEchtrmIneq.push_back(0);
                           }
                       }

                   }
               }

               else
               {
                   if(isVarPrsnt((*itr), var_))
                   {
                       addNumFrEchtrmIneq.push_back(0);
                   }
               }

               maxAddnum.push_back(make_pair(*itr, *max_element(addNumFrEchtrmIneq.begin(),addNumFrEchtrmIneq.end())));
               totalAddnum[*itr] = accumulate(addNumFrEchtrmIneq.begin(), addNumFrEchtrmIneq.end(), 0);
           }


           sort(maxAddnum.begin(), maxAddnum.end(), [=](std::pair<ex, unsigned>& a, std::pair<ex, unsigned>& b){
               return a.second < b.second;}); // sorting in ascending order by maximum number of add in each equation

           unsigned minadd = maxAddnum.begin()->second;

           std::vector<std::pair<ex, unsigned>> temMinAddnum;

           for(auto itr = maxAddnum.begin(); itr != maxAddnum.end(); itr++)
           {
               if((itr->second)==minadd)
                   temMinAddnum.push_back(make_pair(itr->first,totalAddnum[itr->first]));
           }
           if(temMinAddnum.size()>1)
           {
               sort(temMinAddnum.begin(),temMinAddnum.end(),[=](std::pair<ex, unsigned>& a, std::pair<ex, unsigned>& b){
                   return a.second < b.second;}); // sorting in ascending order by total number of add in each equation

               minadd = temMinAddnum.begin()->second;

               std::vector<std::pair<ex, unsigned>> temMinVarnum;

               for(auto itr = temMinAddnum.begin(); itr != temMinAddnum.end(); itr++)
               {
                   if((itr->second)==minadd)
                       temMinVarnum.push_back(make_pair(itr->first,totalAddnum[itr->first]));
               }
               if(temMinVarnum.size()>1)
               {
                   std::vector<std::pair<ex, unsigned>> eq_varnum_pairs;
                   std::map<ex, lst, ex_is_less> eq_var_pairs;

                   for (auto itr = temMinVarnum.begin(); itr != temMinVarnum.end(); ++itr)
                   {
                       for( auto j =  var_.begin(); j !=  var_.end(); j++)
                       {
                           if((itr->first).has(*j))
                               var_in_eq.append(*j);
                       }

                       eq_varnum_pairs.push_back(make_pair(itr->first, nops(var_in_eq)));
                       var_in_eq.remove_all();
                   }
                   sort(eq_varnum_pairs.begin(), eq_varnum_pairs.end(), [=](std::pair<ex, unsigned>& a, std::pair<ex, unsigned>& b){
                       return a.second < b.second;}); // sorting in ascending order by number of var
                   temMinAddnum = eq_varnum_pairs;
               }


           }

           low_var_eq = temMinAddnum.begin()->first;

           if(!nops(var_in_eq))
           {
               for( auto j =  var_.begin(); j !=  var_.end(); j++)
               {
                   if((low_var_eq).has(*j))
                       var_in_eq.append(*j);
               }
           }


           low_var = var_in_eq;


           if(nops(low_var) == 0 || low_var_eq == _ex0)
                return _ex_1;

           exset SysEquClt2;

           unsigned var_soluClt;

           const ex low_var_eq2 = Simplify(collect_common_factors(Factor(low_var_eq)));
           map<unsigned, ex> eqDivider;
           unsigned eqDividerSz = 0;

           if(is_a<mul>(low_var_eq2)) // each factor of Low variable equation is separated by
           {                          // storing into eqDivider. eqDividerSz traces each separated equations.
                ex tem1 =_ex1;
                for(unsigned j = 0; j < nops(low_var_eq2); j++)
                {
                    if(isVarPrsnt((low_var_eq2.op(j)), var_))
                    {
                        if((is_a<power>(low_var_eq2.op(j))&&is_a<add>(low_var_eq2.op(j).op(0)))
                            ||(is_a<add>(low_var_eq2.op(j))))
                        {eqDivider[eqDividerSz] = low_var_eq2.op(j); eqDividerSz++;}
                        else //if(!is_a<numeric>(low_var_eq2.op(j)))
                            tem1 = tem1*low_var_eq2.op(j);
                    }

                }
                if(tem1!=_ex1)
                    {eqDivider[eqDividerSz] = tem1;}
           }
           else if(isVarPrsnt((low_var_eq2), var_))
              eqDivider[0] = low_var_eq2;

           lst low_var_sep;


           eqDividerSz = 0;
           low_var_sep = varOrderByDeg(low_var, eqDivider, eqDividerSz); // solution variables of separated equations are sorted
           lst::const_iterator solu_var = low_var_sep.begin();


           do
           {
                int solu_ext = 0;
                lst tem_one_solu;

                solu_ext = this->one_eq_solutions(Simplify(Numer(eqDivider[eqDividerSz])), (*solu_var)); // attempt 1 for getting solutions

                if(nops(one_eq_solu))
                {
                    lst one_eq_solu_test = {};
                    for(auto itr1 = one_eq_solu.begin(); itr1 != one_eq_solu.end(); itr1++)
                    {
                        if(simplifyRecur((Numer(eqDivider[eqDividerSz])).subs(*itr1))==_ex0)
                            one_eq_solu_test.append(*itr1);
                    }

                    one_eq_solu = one_eq_solu_test;
                    if(nops(one_eq_solu))
                        solu_ext =1;
                    else
                    {
                        solu_ext = 0;
                    }
                }
                if(solu_ext) // avoiding pole error
                {
                    lst one_solu_fltr={};
                    for(auto itr1 = one_eq_solu.begin(); itr1 != one_eq_solu.end(); itr1++)
                    {
                        isSolu = true;

                        for(auto j = sysequ.begin(); j != sysequ.end();)
                        {
                            denoma = Denom(*j);
                            if(denoma != _ex1 && denoma.has((*solu_var)))
                            {
                                try
                                {
                                    subs_simp = (simplify((denoma).subs(*itr1)));
                                    if((is_a<numeric>(subs_simp)  || (is_a<numeric>(subs_simp) && subs_simp.info(info_flags::positive) && subs_simp < Gtolerance) ||
                                        (is_a<numeric>(subs_simp) && subs_simp.info(info_flags::negative) && -subs_simp < Gtolerance)))
                                    {
                                        isSolu = false;
                                        j = sysequ.end();
                                    }
                                    else
                                        j++;

                                }
                                catch(GiNaC::pole_error){isSolu = false; j = sysequ.end();}
                            }
                            else
                                j++;

                        }

                        if(isSolu)
                        {
                            one_solu_fltr.append(*itr1);
                        }
                    }

                    if(nops(one_solu_fltr)==0)
                    {
                        one_eq_solu.remove_all();
                        solu_ext = 0;
                    }
                    else
                        one_eq_solu = one_solu_fltr;

                }

               if(nops(one_eq_solu) != 0)
               {

                  if(sysequ.find(low_var_eq)!=sysequ.end())
                    sysequ.erase(low_var_eq);
                   if(nops(one_eq_solu) != 0)
                     tem_one_solu = one_eq_solu;

                   one_eq_solu.remove_all();

                   var_soluClt = nops(soluClt);

                   for(lst::const_iterator k = tem_one_solu.begin(); k != tem_one_solu.end(); k++) // searching final solutions taking each solution of low variable equation
                   {

                       soluClt.append((*k).lhs() == Simplify(Factor(Numer((*k).rhs()), factor_all))/Simplify(Factor(Denom((*k).rhs()), factor_all)));

                       for(auto j = sysequ.begin(); j != sysequ.end(); j++)
                       {
                            try
                            {
                                subs_simp = Numer(simplify(Numer(simplify(*j)).subs(simplify(*k)))); // GiNaC Bugs: wrong numer
                            }                                                                        // numer(1/2*b/((1 + 3*sqrt(a)*(-2 + 9*a)^(-1/2))*(-1 + 3*sqrt(a)*(-2 + 9*a)^(-1/2))));
                            catch(GiNaC::pole_error){return _ex_1;}                                  //  -> -(-2+9*a)*b
                            catch(cln::runtime_exception){return _ex_1;}                             // To avoid Bugs: get numer after simplify

                            if( !is_a<numeric>(subs_simp) ||
                               (is_a<numeric>(subs_simp) && subs_simp.info(info_flags::positive) && subs_simp > Gtolerance) ||
                               (is_a<numeric>(subs_simp) && subs_simp.info(info_flags::negative) && -subs_simp > Gtolerance))
                                SysEquClt2.insert(subs_simp);

                       }

                       SysEquCltEnv = SysEquClt2;

                       const ex ret = this->sysequ_red(SysEquClt2,  var_);

                       if(nops(soluClt) - (var_soluClt) == 1)
                       {
                            soluClt.remove_last();
                       }

                        SysEquClt2.clear();
                   }
                }

                solu_var++;

                if(eqDivider.size() > (eqDividerSz + 1) && solu_var == low_var_sep.end()) // the solution process is restarted
                {                                                                         // for next equation in eqDivider.
                    eqDividerSz++;

                    low_var_sep = varOrderByDeg(low_var, eqDivider, eqDividerSz);
                    if(nops(low_var_sep != 0))
                        solu_var = low_var_sep.begin();
                    else
                        solu_var = low_var_sep.end();

                }
                else if(eqDivider.size() > (eqDividerSz + 1) && SysEquCltEnv.empty())  // the solution process for current equation is stopped and
                {                                                                      // we do not get the solutions for other variables.
                    eqDividerSz++;                                                    // the solution process is restarted
                    low_var_sep = varOrderByDeg(low_var, eqDivider, eqDividerSz);      // for next equation in eqDivider.
                    if(nops(low_var_sep != 0))
                        solu_var = low_var_sep.begin();
                    else
                        solu_var = low_var_sep.end();
                }
                else if(SysEquCltEnv.empty() && ((is_a<power>(eqDivider[eqDividerSz])&&is_a<add>(eqDivider[eqDividerSz].op(0)))
                                                ||(is_a<add>(eqDivider[eqDividerSz])))) // reaches at maximum equation and the solution
                    {solu_var = low_var_sep.end();}                                 // process for current equation is stopped.

            }while(solu_var != low_var_sep.end());
        }

        catch(GiNaC::pole_error){return _ex1;}
        catch (cln::runtime_exception){return _ex1; }

        catch(std::invalid_argument){return _ex1; }
        catch(std::out_of_range){return _ex1;}



    }

    return _ex1;
}

/** All solutions take the shorten forms by repeated substitutions.  **/
exsetlst solvec::solu_subs(exsetlst solu_)
{
    exsetlst solu_wt_subs;
    lst tem;
    ex xprev, xnow;
    bool den_zero;

   /// solutions are substituted by other solutions
    for(exsetlst::const_iterator it = solu_.begin(); it != solu_.end(); it++)
    {
        unsigned i = 0;
        den_zero = false;
        exset lst_exset;
        lst exset_lst;
        do
        {

            xnow = ((*it).op(i)).rhs();

            do
            {
                xprev = (xnow);

                try
                {
                    xnow = simplify((xprev)).subs((*it), subs_options::algebraic);
                }
                catch(GiNaC::pole_error){den_zero = true;}

            } while(xprev != xnow && den_zero == false);
            if(den_zero == false)
                tem.append(((*it).op(i)).lhs() == xnow);
         i = i + 1;
        }while(den_zero == false && i < nops(*it));

        if(den_zero == false)
        {
            lst_exset.insert(tem.begin(), tem.end());
            for(exset::const_iterator it = lst_exset.begin(); it != lst_exset.end(); it++)
            {
                exset_lst.append(*it);
            }
            solu_wt_subs.insert(exset_lst);
            exset_lst.remove_all();
            lst_exset.clear();
        }
        tem.remove_all();
    }

    /// the presence of sub-solutions in each solutions are checked. If sub-solution present, the solutions
    /// are removed.
    vector<std::pair<lst, unsigned>> eq_nops_pairs;
    for(exsetlst::const_iterator it = solu_wt_subs.begin(); it != solu_wt_subs.end(); it++)
    {
       eq_nops_pairs.push_back(make_pair(*it, nops(*it)));
    }
    sort(eq_nops_pairs.begin(), eq_nops_pairs.end(), [=](std::pair<lst, unsigned>& a, std::pair<lst, unsigned>& b){
    return a.second < b.second;});

    for(vector<std::pair<lst, unsigned>>::const_iterator it1 = eq_nops_pairs.begin(); it1 != (eq_nops_pairs.end()); it1++)
    {
        vector<std::pair<lst, unsigned>>::const_iterator it3 = it1 + 1;
        for(auto it4 = (it3); it4 != (eq_nops_pairs).end(); it4++)
        {
            unsigned lst_num = 0;
            for(lst::const_iterator it2 = (it1->first).begin(); it2 != (it1->first).end(); it2++)
            {
                for(lst::const_iterator it5 = ((it4)->first).begin(); it5 != (it4->first).end(); it5++)
                {
                    if((*it2) == (*it5))
                        lst_num = lst_num + 1;
                }
            }

            if(lst_num == nops(it1->first))
            {
                solu_wt_subs.erase((it4->first));
            }

        }

    }

    return solu_wt_subs;
}
/////////////////////////////////////////////////////////////

exsetlst solvec::operator()(const lst & equ_, const lst& var_)
{
    ex numera, denoma, subs_imp;
    exset _var, sysequ;
    exset equ_set;
    exsetlst temsolu_setlst,solu_setlst;
    map<unsigned, exset> eqDivider, temeqDivider;

    try
    {

        for (lst::const_iterator i = var_.begin(); i != var_.end(); ++i)
        {
           _var.insert(*i);
        }

        for (lst::const_iterator i = equ_.begin(); i != equ_.end(); ++i)
        {
            numera = simplify(Numer(*i));

            for( auto itr = _var.begin(); itr != _var.end(); itr++ )
            {
                funcpresent funcpresent(*itr);
                funcpresent(numera);
                if(funcpresent.funcpresence || funcpresent.varInPow) // presence of any functions with solution variables,
                    return solu_setlst;                              // presence of solution variables in power are not supported.
            }

           if(numera != _ex0)
           {
                equ_set.insert(numera);
           }

        }

        if(!equ_set.empty())
            this->sysequ_red(equ_set, _var);

        if(!solu.empty())
            temsolu_setlst = this->solu_subs(solu);

        bool notSolu;
        for(auto itr = temsolu_setlst.begin();itr!=temsolu_setlst.end();itr++)
        {
            notSolu = false;
            for (auto itr1 = equ_set.begin();itr1 != equ_set.end();)
            {
                if(simplifyRecur((*itr1).subs(*itr))!=_ex0)
                {
                    notSolu =true;
                    itr1 = equ_set.end();
                }
                else
                {
                    itr1++;
                }
            }

            if(!notSolu)
                solu_setlst.insert(*itr);
        }

        return solu_setlst;
    }

    catch(GiNaC::pole_error){return solu_setlst;}
    catch (cln::runtime_exception){return solu_setlst;}

    catch(std::invalid_argument){return solu_setlst;}
    catch(std::out_of_range){return solu_setlst;}
    //catch(std::runtime_error){return solu_setlst;}
    //catch(std::range_error){return solu_setlst;}
    //catch(std::logic_error){return solu_setlst;}
    //catch(std::domain_error){return solu_setlst;}

}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////







