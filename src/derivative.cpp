
/** @file derivative.cpp
 *
 *  Implementation of GiNaCDE's derivative. */



#include <ginac/ginac.h>
#include "utility.h"
#include "derivative.h"
#include "integrate.h"
#include "simplify.h"
#include "inifcns.h"

using namespace std;
using namespace GiNaC;

/** removes zero time Diff, Integrate terms **/
ex zero_order_rem(const ex& expr_)
{
    ex _y=expr_;
    exmap expr;
    expr[Integrate(wild(0), wild(1), 0)] = wild(0);
    expr[Diff(wild(0), wild(1), 0)] = wild(0);
    ex xprev;
    do
    {
        xprev = _y;
        _y = _y.subs(expr, subs_options::algebraic);
    } while(xprev != _y);
    return _y;
}
//////////////////////////////////////////////////////////////////////////////////


ex subs_IntDiffargu::operator()(const ex& _e)
{
    if (is_ex_the_function(_e, Integrate)
        || (is_ex_the_function(_e, Diff)))
    {
        str = "IntDiffargu_" + to_string(i);

        _expr = reader(str);
        _temexpr = _e.subs(_e.op(0) == _expr, subs_options::algebraic);

        if (_temexpr.has(_var))
        {
            IntDiffargu_map[_expr] = _e.op(0);
            i = i + 1;
            return _temexpr;
        }
        else
        {
            symclt = symbols(_e.op(0));
            for (auto it1 = symclt.begin(); it1 != symclt.end();)
            {

                auto it2 = var_depend.find(*it1);
                if (it2 != var_depend.end())
                {
                    auto it = find(var_depend[*it1].begin(), var_depend[*it1].end(), _var);
                    if (it != var_depend[*it1].end())
                    {
                        var_depend[_expr] = {_var};
                        it1 = symclt.end();
                    }
                    else
                    {
                         it1++;
                    }
                }
                else
                {
                     it1++;
                }
            }

            IntDiffargu_map[_expr] = _e;
            i = i + 1;
            return _expr;
        }
    }

    else
    {
        //doing derivative of conjugate function

        if( nops(_e) == 1 && _e.op(0) == conjugate(_e) )
        {
            ex _tem = _e;
            str = "Conjuargu_" + to_string(j);

            _expr = reader(str);
            _temexpr = _e.subs(( _tem ) == _expr, subs_options::algebraic);

            symclt = symbols(_tem);
            for (auto it1 = symclt.begin(); it1 != symclt.end(); )
            {
                    if (*it1 == _var)
                    {
                        var_depend[_expr] = {_var};
                        it1 = symclt.end();
                    }
                    else
                    {
                        it1++;
                    }
            }
            for (auto it1 = symclt.begin(); it1 != symclt.end(); )
            {
                auto it2 = var_depend.find(*it1);
                if (it2 != var_depend.end())
                {
                    auto it = find(var_depend[*it1].begin(), var_depend[*it1].end(), _var);
                    if (it != var_depend[*it1].end())
                    {
                        var_depend[_expr] = {_var};
                        it1 = symclt.end();
                    }
                    else
                    {
                        it1++;
                    }
                }
                else
                    it1++;
            }


            Conjurgu_map[_expr] = ( _tem );
            j = j + 1;

            return _temexpr;

        }
        else
        {
            return _e.map(*this);
        }

    }

}
//////////////////////////////////////////////////////////////////////////////////


ex pdiff(const ex& expr_, const ex& var_, const ex& order_)
{
    if(!is_a<symbol>(var_) || !is_a<numeric>(order_) || ex_to<numeric>(order_).is_negative())
        return expr_;
    ex _expr;
    subs_IntDiffargu subs_IntDiffargu(0, 0, var_);

    _expr = subs_IntDiffargu(expr_);
    exmap IntDiffargu = subs_IntDiffargu.IntDiffargu_map;
    exmap Conjuargu = subs_IntDiffargu.Conjurgu_map;
    exset setex;
    exmap mexsymb;
    setex = symbols(_expr);

    for (auto it1 = setex.begin(); it1 != setex.end(); it1++)
    {
        auto it2 = var_depend.find(*it1);
        if (it2 != var_depend.end())
        {
            auto it = find(var_depend[*it1].begin(), var_depend[*it1].end(), var_);
            if (it != var_depend[*it1].end())
            {
                mexsymb[*it1] = Diff(*it1, var_, 0);
            }
        }
    }

    if(!mexsymb.empty())
    {
        _expr = _expr.subs(mexsymb, subs_options::algebraic);
    }

    if(!(order_).info(info_flags::integer))
    {
        dorat dorat; // avoiding double form
        dorat.set();
        if(!dorat(order_).info(info_flags::integer))
            return _FAIL;
        _expr = diff(_expr, ex_to<symbol>(var_), ex_to<numeric>(dorat(order_)).to_int());
    }
    else
        _expr = diff(_expr, ex_to<symbol>(var_), ex_to<numeric>((order_)).to_int());

   if (!IntDiffargu.empty())
   {
       _expr = _expr.subs(IntDiffargu, subs_options::algebraic);
       _expr = zero_order_rem(_expr);
       for (auto it = IntDiffargu.begin(); it != IntDiffargu.end(); it++)
       {
          var_depend.erase(it->first);
       }
   }

   if( !Conjuargu.empty() )
   {

       _expr = _expr.subs(Conjuargu, subs_options::algebraic);
       _expr = zero_order_rem(_expr);
       for (auto it = Conjuargu.begin(); it != Conjuargu.end(); it++)
       {
          var_depend.erase(it->first);
       }
   }


    return zero_order_rem(_expr);
}
//////////////////////////////////////////////////////////////////////////////////

evaluatec evaluateg;
ex evaluatec::operator()(const ex& _e)
{
    if (is_ex_the_function(_e, Diff) && !is_ex_the_function(_e.op(0), Diff))
       return pdiff(_e.op(0), _e.op(1), _e.op(2));

    if (is_ex_the_function(_e, Integrate) && !is_ex_the_function(_e.op(0), Integrate))
    {
        if(_e.op(2) != _ex0)
        {
            _tem = integrate(_e.op(0), _e.op(1));

            int extonum;
            if(!(_e.op(2)).info(info_flags::integer))
            {
                dorat dorat; // avoiding double form
                dorat.set();
                if(!dorat(_e.op(2)).info(info_flags::integer))
                    return _FAIL;
                extonum=ex_to<numeric>(dorat(_e.op(0))).to_int();
            }
            else
                extonum=ex_to<numeric>(_e.op(0)).to_int();

            for(int i = 1; i < extonum; i++)
                _tem = integrate(_tem, _e.op(1));

            return _tem;
        }
        else
            return _e.op(0);

    }

    return _e.map(*this);
}
//////////////////////////////////////////////////////////////////////////////////

ex evaluate(const ex& expr_)
{
    ex _y=Simplify(expand(expr_));
    //ex _y=(expand(expr_));

    ex _xprev;
    do
    {
        _xprev = _y;
        _y = evaluateg(_y);
    } while(_xprev != _y);

    //return Simplify(expand(_y));
    return _y;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


