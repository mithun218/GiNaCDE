
/** @file inifcns.cpp
 *
 *  Implementation of some usefull funtions (acot, asec, acsc, cot, sec, csc,
    acoth, asech, acsch, coth, sech, csch, sn, cn, dn, Diff, Integrate) */



#include <ginac/ginac.h>
#include "inifcns.h"
#include "utility.h"

using namespace std;
using namespace GiNaC;



const numeric *_num_1_p = (const numeric *)&dynallocate<numeric>(-1);
const numeric *_num0_p = (const numeric *)&dynallocate<numeric>(0);
const numeric *_num1_p = (const numeric *)&dynallocate<numeric>(1);
 //-------------------------------------------------



 //////////
// inverse cotangent (arc cotangent)
//////////

static ex acot_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return atan(ex_to<numeric>(x).inverse());

        return acot(x).hold();
}

static ex acot_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {

                if (x.is_zero())
                        return _ex1_2*Pi;

                if (x.is_equal(_ex1))
                        return _ex1_4*Pi;

                if (x.is_equal(_ex_1))
                        return _ex_1_4*Pi;

                if (x.is_equal(I) || x.is_equal(-I))
                        throw (pole_error("acot_eval(): logarithmic pole",0));

        if (!x.info(info_flags::crational))
            return atan(ex_to<numeric>(x).inverse());

                if (x.info(info_flags::negative))
                        return -acot(-x);
        }

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);
                if (is_ex_the_function(x, cot))
                        return t;
        }
        return acot(x).hold();
}

static ex acot_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);
        return mul(power(_ex1+power(x,_ex2), _ex_1), _ex_1);
}

static ex acot_series(const ex &arg,
                      const relational &rel,
                      int order,
                      unsigned options)
{
        return _ex1_2*Pi -(atan(arg)).series(rel, order, options);
}

REGISTER_FUNCTION(acot, eval_func(acot_eval).
                        evalf_func(acot_evalf).
                        derivative_func(acot_deriv).
                        series_func(acot_series).
                        set_name("arccot", "\\operatorname{arccot}"));


//////////
// inverse secant (arc secant)
//////////

static ex asec_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return acos(ex_to<numeric>(x).inverse());

        return asec(x).hold();
}

static ex asec_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                //numeric num = ex_to<numeric>(x);
                if (x.is_zero())
                        throw (pole_error("asec_eval(): logarithmic pole",0));
                if (x.is_equal(_ex1))
                        return _ex0;
                if (x.is_equal(_ex_1))
                        return Pi;
                if (!x.info(info_flags::crational))
                        return acos(ex_to<numeric>(x).inverse());
        }


        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);
                if (is_ex_the_function(x, sec))
                        return t;
        }
        return asec(x).hold();
}

static ex asec_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);
        return power(mul(x, power(_ex_1 + power(x,_ex2), _ex1_2)), _ex_1);
}

REGISTER_FUNCTION(asec, eval_func(asec_eval).
                        evalf_func(asec_evalf).
                        derivative_func(asec_deriv).
                        set_name("arcsec", "\\operatorname{arcsec}"));


//////////
// inverse cosecant (arc cosecant)
//////////

static ex acsc_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return asin(ex_to<numeric>(x).inverse());

        return acsc(x).hold();
}

static ex acsc_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                if (x.is_zero())
                        throw (pole_error("acsc_eval(): logarithmic pole",0));
                if (x.is_equal(_ex1))
                        return Pi/_ex2;
                if (x.is_equal(_ex_1))
                        return -Pi/_ex2;
                if (!x.info(info_flags::crational))
                        return asin(ex_to<numeric>(x).inverse());
        }

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);
                if (is_ex_the_function(x, csc))
                        return t;
        }
        return acsc(x).hold();
}

static ex acsc_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);
        return -power(mul(x, power(_ex_1 + power(x,_ex2), _ex1_2)), _ex_1);
}

REGISTER_FUNCTION(acsc, eval_func(acsc_eval).
                        evalf_func(acsc_evalf).
                        derivative_func(acsc_deriv).
                        set_name("arccsc", "\\operatorname{arccsc}"));



//////////
// cotangent (trigonometric function)
//////////

static ex cot_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                if (ex_to<numeric>(x).is_zero())
                    throw (pole_error("cot_evalf(): simple pole",1));
                return tan(ex_to<numeric>(x)).inverse();
        }

        return cot(x).hold();
}

static ex cot_eval(const ex & x)
{
        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                 //cot(acot(x)) -> x
                if (is_ex_the_function(x, acot))
                        return t;

                // cot(asin(x)) -> sqrt(1-x^2)/x
                if (is_ex_the_function(x, asin))
                        return sqrt(_ex1-power(t,_ex2))/t;

                // cot(acos(x)) -> x/sqrt(1-x^2)
                if (is_ex_the_function(x, acos))
                        return t*power(_ex1-power(t,_ex2),_ex_1_2);

                // cot(atan(x)) -> 1/x;
                if (is_ex_the_function(x, atan))
                        return _ex1/t;

                // cot(acsc(x)) -> sqrt(x^2-1)
                if (is_ex_the_function(x, acsc))
                        return sqrt(power(t,_ex2)-_ex1);

                // cot(asec(x)) -> 1/sqrt(x^2-1)
                if (is_ex_the_function(x, asec))
                        return power(power(t,_ex2)-_ex1,_ex_1_2);

                // cot(atan2(y,x)) -> x/y;
                if (is_ex_the_function(x, atan2)) {
                        const ex &t1 = x.op(1);
                        return t1/t;
                }
        }

        // cot(float) -> float
        if (x.info(info_flags::numeric) && !x.info(info_flags::crational)) {
        if (ex_to<numeric>(x).is_zero())
            throw (pole_error("cot_eval(): simple pole",1));
        return tan(ex_to<numeric>(x)).inverse();
        }
    if (!(((ex)x/Pi*2).info(info_flags::odd)))
    {
        ex res = tan(x);
        if (not is_ex_the_function(res, tan) && not is_ex_the_function(_ex_1*res, tan)) {
            if (not res.is_zero()) {
                return tan((ex)(Pi/2-x));
                }
            else
                throw (pole_error("cot_eval(): simple pole",1));
        }
    }
        // Reflection at Pi/2
        const ex ExOverPi = x/Pi;
        if(is_exactly_a<numeric>(ExOverPi))
        {
            ex coef_pi = x.coeff(Pi).expand();
            if (is_exactly_a<numeric>(coef_pi))
            {
                const numeric c = ex_to<numeric>(coef_pi);
                if (c.is_rational())
                {
                    const numeric num = c.numer();
                    const numeric den = c.denom();
                    const numeric rm = mod(num,den);
                    if (rm.mul(2) == den)
                    {
                        return _ex0;
                    }
                }
            }
        }
        return cot(x).hold();
}

static ex cot_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx cot(x) -> -1-cot(x)^2;
        return (_ex_1-power(cot(x),_ex2));
}

// Ref.: http://dlmf.nist.gov/4.21.E40
static ex cot_real_part(const ex & x)
{
        ex a = GiNaC::real_part(mul(x, _ex2));
        ex b = GiNaC::imag_part(mul(x, _ex2));
        return sin(a)/(cosh(b) - cos(a));
}

static ex cot_imag_part(const ex & x)
{
        ex a = GiNaC::real_part(mul(x, _ex2));
        ex b = GiNaC::imag_part(mul(x, _ex2));
        return sinh(b)/(cosh(b) - cos(a));
}

static ex cot_series(const ex &x,
                     const relational &rel,
                     int order,
                     unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        // method:
        // Taylor series where there is no pole falls back to tan_deriv.
        // On a pole simply expand cos(x)/sin(x).
        const ex x_pt = x.subs(rel, subs_options::no_pattern);
        if (!(2*x_pt/Pi).info(info_flags::even))
                throw do_taylor();  // caught by function::series()
        // if we got here we have to care for a simple pole
        return (cos(x)/sin(x)).series(rel, order, options);
}

static ex cot_conjugate(const ex & x)
{
        // conjugate(tan(x))==1/tan(conjugate(x))
        return power(tan(x.conjugate()), _ex_1);
}

REGISTER_FUNCTION(cot, eval_func(cot_eval).
                       evalf_func(cot_evalf).
                       derivative_func(cot_deriv).
                       series_func(cot_series).
                       real_part_func(cot_real_part).
                       imag_part_func(cot_imag_part).
                       conjugate_func(cot_conjugate).
                       latex_name("\\cot"));

//////////
// secant (trigonometric function)
//////////

static ex sec_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return cos(ex_to<numeric>(x)).inverse();

        return sec(x).hold();
}

static ex sec_eval(const ex & x)
{
        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // sec(asec(x)) -> x
                if (is_ex_the_function(x, asec))
                        return t;

                // sec(asin(x)) -> 1/sqrt(1-x^2)
                if (is_ex_the_function(x, asin))
                        return power(_ex1-power(t,_ex2),_ex_1_2);

                // sec(acos(x)) -> 1/x
                if (is_ex_the_function(x, acos))
                        return _ex1/t;

                // sec(atan(x)) -> sqrt(x^2+1)
                if (is_ex_the_function(x, atan))
                        return sqrt(power(t,_ex2)+_ex1);

                // sec(acot(x)) -> sqrt(x^2+1)/x
                if (is_ex_the_function(x, acot))
                        return sqrt(power(t,_ex2)+_ex1)/t;

                // sec(acsch(x)) -> x/sqrt(x^2-1)
                if (is_ex_the_function(x, acsc))
                        return t*power(power(t,_ex2)-_ex1,_ex_1_2);

                // sec(atan2(y,x)) -> sqrt(x^2+y^2)/x
                if (is_ex_the_function(x, atan2)) {
                        const ex &t1 = x.op(1);
                        return sqrt(power(t1,_ex1)+power(t,_ex2))/t1;
                }

        }

        // sec() is even
        if (x.info(info_flags::negative))
                return sec(-x);

        // sec(float) -> float
        if (is_exactly_a<numeric>(x) && !x.info(info_flags::crational)) {
                return cos(ex_to<numeric>(x)).inverse();
        }

        // Handle simplification via sin
        ex res = cos(x);
        if (not is_ex_the_function(res, cos) && not is_ex_the_function(_ex_1*res, cos)) {
                if (res.is_zero())
                        throw (pole_error("sec_eval(): simple pole",1));
                else
                        return power(res, _ex_1);
        }

        // cos has reflected also the argument so take it
        if (is_ex_the_function(res, cos))
                return sec(res.op(0)).hold();
        else
                return -sec((-res).op(0)).hold();
}

static ex sec_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx sec(x) -> sec(x)*tan(x);
        return sec(x)*tan(x);
}

static ex sec_real_part(const ex & x)
{
        ex a = GiNaC::real_part(x);
        ex b = GiNaC::imag_part(x);
        return cos(a)*cosh(b)/(sin(a)*sin(a)*sinh(b)*sinh(b) \
            + cos(a)*cos(a)*cosh(b)*cosh(b));
}

static ex sec_imag_part(const ex & x)
{
        ex a = GiNaC::real_part(x);
        ex b = GiNaC::imag_part(x);
        return sin(a)*sinh(b)/(sin(a)*sin(a)*sinh(b)*sinh(b) \
            + cos(a)*cos(a)*cosh(b)*cosh(b));
}

static ex sec_series(const ex &x,
                     const relational &rel,
                     int order,
                     unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        return (_ex1/cos(x)).series(rel, order, options);
}

static ex sec_conjugate(const ex & x)
{
        // conjugate(tan(x))==1/tan(conjugate(x))
        return power(cos(x.conjugate()), _ex_1);
}

REGISTER_FUNCTION(sec, eval_func(sec_eval).
                       evalf_func(sec_evalf).
                       derivative_func(sec_deriv).
                       series_func(sec_series).
                       real_part_func(sec_real_part).
                       imag_part_func(sec_imag_part).
                       conjugate_func(sec_conjugate).
                       latex_name("\\sec"));

//////////
// cosecant (trigonometric function)
//////////

static ex csc_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                if (ex_to<numeric>(x).is_zero())
                        throw (pole_error("csc_eval(): simple pole",1));
                return sin(ex_to<numeric>(x)).inverse();
        }

        return csc(x).hold();
}

static ex csc_eval(const ex & x)
{

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // csc(acsc(x)) -> x
                if (is_ex_the_function(x, acsc))
                        return t;

                // csc(asin(x)) -> 1/x
                if (is_ex_the_function(x, asin))
                        return _ex1/t;

                // csc(acos(x)) -> 1/sqrt(1-x^2)
                if (is_ex_the_function(x, acos))
                        return power(_ex1-power(t,_ex2),_ex_1_2);

                // csc(atan(x)) -> sqrt(x^2+1)/x;
                if (is_ex_the_function(x, atan))
                        return sqrt(power(t,_ex2)+_ex1)/t;

                // csc(acot(x)) -> sqrt(x^2+1);
                if (is_ex_the_function(x, acot))
                        return sqrt(power(t,_ex2)+_ex1);

                // csc(asec(x)) -> x/sqrt(x^2-1)
                if (is_ex_the_function(x, asec))
                        return t*power(power(t,_ex2)-_ex1,_ex_1_2);

                // csc(atan2(y,x)) -> sqrt(x^2+y^2)/y
                if (is_ex_the_function(x, atan2)) {
                        const ex &t1 = x.op(1);
                        return sqrt(power(t1,_ex1)+power(t,_ex2))/t;
                }

        }

        // csc(float) -> float
        if (is_exactly_a<numeric>(x) && !x.info(info_flags::crational)) {
                if (ex_to<numeric>(x).is_zero())
                        throw (pole_error("csc_eval(): simple pole",1));
                return sin(ex_to<numeric>(x)).inverse();
        }

        // Handle simplification via sin
        ex res = sin(x);
        if (not is_ex_the_function(res, sin) && not is_ex_the_function(_ex_1*res, sin)) {
                if (res.is_zero())
                        throw (pole_error("csc_eval(): simple pole",1));
                else
                        return power(res, _ex_1);
        }
        // sin has reflected also the argument so take it
        if (is_ex_the_function(res, sin))
                return csc(res.op(0)).hold();
        else
                return -csc((-res).op(0)).hold();
}

static ex csc_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx cot(x) -> -1-cot(x)^2;
        return -csc(x)*cot(x);
}

static ex csc_real_part(const ex & x)
{
        ex a = GiNaC::real_part(x);
        ex b = GiNaC::imag_part(x);
        return sin(a)*cosh(b)/(sin(a)*sin(a)*cosh(b)*cosh(b) \
            + cos(a)*cos(a)*sinh(b)*sinh(b));
}

static ex csc_imag_part(const ex & x)
{
        ex a = GiNaC::real_part(x);
        ex b = GiNaC::imag_part(x);
        return -cos(a)*sinh(b)/(sin(a)*sin(a)*cosh(b)*cosh(b) \
            + cos(a)*cos(a)*sinh(b)*sinh(b));
}

static ex csc_series(const ex &x,
                     const relational &rel,
                     int order,
                     unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        return (_ex1/sin(x)).series(rel, order, options);
}

static ex csc_conjugate(const ex & x)
{
        // conjugate(tan(x))==1/tan(conjugate(x))
        return power(sin(x.conjugate()), _ex_1);
}

REGISTER_FUNCTION(csc, eval_func(csc_eval).
                       evalf_func(csc_evalf).
                       derivative_func(csc_deriv).
                       series_func(csc_series).
                       real_part_func(csc_real_part).
                       imag_part_func(csc_imag_part).
                       conjugate_func(csc_conjugate).
                       latex_name("\\csc"));

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////
// hyperbolic sine (trigonometric function)
//////////



//////////
// hyperbolic cotangent (trigonometric function)
//////////

static ex coth_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return tanh(ex_to<numeric>(x)).inverse();

        return tanh(x).hold();
}

static ex coth_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {

                // coth(0) -> zoo
                if (x.is_zero())
                        throw pole_error("coth_eval(): simple pole",1);

                // coth(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return tanh(ex_to<numeric>(x)).inverse();

                // coth() is odd
                if (x.info(info_flags::negative))
                        return -coth(-x);
        }


        ex xoverpi = x/Pi;
        if (is_exactly_a<numeric>(xoverpi)) {
                numeric nxopi = ex_to<numeric>(xoverpi);
                if (nxopi.real().is_zero()) { // coth(I*x) -> I*cot(x);
                        numeric xoverpiI = (nxopi*2)/I;
                        if (not xoverpiI.is_integer())
                                return -I*cot(x/I);
                        else if (xoverpiI.is_odd())
                                return _ex0;
                        else
                                throw pole_error("coth_eval(): simple pole",1);
                }
                else
                        return coth(x).hold();
        }

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // coth(acoth(x)) -> x
                if (is_ex_the_function(x, acoth))
                        return t;

                // coth(asinh(x)) -> sqrt(1+x^2)/x
                if (is_ex_the_function(x, asinh))
                        return power(_ex1+power(t,_ex2),_ex1_2)/t;

                // coth(acosh(x)) -> x/(sqrt(x-1)*sqrt(x+1))
                if (is_ex_the_function(x, acosh))
                        return t/sqrt(t-_ex1)/sqrt(t+_ex1);
        }

        return coth(x).hold();
}

static ex coth_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx coth(x) -> 1-coth(x)^2
        return 1-power(coth(x),_ex2);
}

static ex coth_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        // method:
        // Taylor series where there is no pole falls back to tanh_deriv.
        // On a pole simply expand sinh(x)/cosh(x).
        const ex x_pt = x.subs(rel, subs_options::no_pattern);
        if (!(2*I*x_pt/Pi).info(info_flags::even))
                throw do_taylor();  // caught by function::series()
        // if we got here we have to care for a simple pole
        return (cosh(x)/sinh(x)).series(rel, order, options);
}

static ex coth_real_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return mul(sinh(a), cosh(a)) / (power(sin(b), _ex2) + power(sinh(a), _ex2));
}

static ex coth_imag_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return -mul(sin(b), cos(b)) / (power(sin(b), _ex2) + power(sinh(a), _ex2));
}

static ex coth_conjugate(const ex & x)
{
        return coth(x.conjugate());
}

REGISTER_FUNCTION(coth, eval_func(coth_eval).
                        evalf_func(coth_evalf).
                        derivative_func(coth_deriv).
                        series_func(coth_series).
                        real_part_func(coth_real_part).
                        imag_part_func(coth_imag_part).
                        conjugate_func(coth_conjugate).
                        latex_name("\\coth"));

//////////
// hyperbolic secant (trigonometric function)
//////////

static ex sech_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return cosh(ex_to<numeric>(x)).inverse();

        return sech(x).hold();
}

static ex sech_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {

                // sech(0) -> 1
                if (x.is_zero())
                        return _ex1;

                // sech(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return cosh(ex_to<numeric>(x)).inverse();

                // sech() is even
                if (x.info(info_flags::negative))
                        return sech(-x);
        }

        ex xoverpi = x/Pi;
        if (is_exactly_a<numeric>(xoverpi) &&
                ex_to<numeric>(xoverpi).real().is_zero())
                return sec(x/I);

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // sech(asech(x)) -> x
                if (is_ex_the_function(x, asech))
                        return t;

                // sech(asinh(x)) -> 1/sqrt(1+x^2)
                if (is_ex_the_function(x, asinh))
                        return power(_ex1+power(t,_ex2),_ex_1_2);

                // sech(acosh(x)) -> 1/x
                if (is_ex_the_function(x, acosh))
                        return power(t, _ex_1);
        }

        return sech(x).hold();
}

static ex sech_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx sech(x) -> -sech(x)*tanh(x)
        return -mul(sech(x), tanh(x));
}

static ex sech_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        // method:
        // Taylor series where there is no pole falls back to tanh_deriv.
        // On a pole simply expand sinh(x)/cosh(x).
        const ex x_pt = x.subs(rel, subs_options::no_pattern);
        if (!(2*I*x_pt/Pi).info(info_flags::odd))
                throw do_taylor();  // caught by function::series()
        // if we got here we have to care for a simple pole
        return (_ex1/cosh(x)).series(rel, order, options);
}

static ex sech_real_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return mul(cos(b), cosh(a)) / (power(mul(sin(b), sinh(a)), _ex2) + power(mul(cos(b), cosh(a)), _ex2));
}

static ex sech_imag_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return -mul(sin(b), sinh(a)) / (power(mul(sin(b), sinh(a)), _ex2) + power(mul(cos(b), cosh(a)), _ex2));
}

static ex sech_conjugate(const ex & x)
{
        return sech(x.conjugate());
}

REGISTER_FUNCTION(sech, eval_func(sech_eval).
                        evalf_func(sech_evalf).
                        derivative_func(sech_deriv).
                        series_func(sech_series).
                        real_part_func(sech_real_part).
                        imag_part_func(sech_imag_part).
                        conjugate_func(sech_conjugate).
                        latex_name("\\operatorname{sech}"));

//////////
// hyperbolic secant (trigonometric function)
//////////

static ex csch_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return sinh(ex_to<numeric>(x)).inverse();

        return csch(x).hold();
}

static ex csch_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {

                // csch(0) -> zoo
                if (x.is_zero())
                        throw pole_error("csch_eval(): simple pole",1);

                // csch(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return sinh(ex_to<numeric>(x)).inverse();

                // csch() is odd
                if (x.info(info_flags::negative))
                        return -csch(-x);
        }

        ex xoverpi = x/Pi;
        if (is_exactly_a<numeric>(xoverpi) &&
                ex_to<numeric>(xoverpi).real().is_zero())
                return -I*csc(x/I);


        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // sech(asech(x)) -> x
                if (is_ex_the_function(x, acsch))
                        return t;

                // sech(acosh(x)) -> 1/sqrt(x^2-1)
                if (is_ex_the_function(x, asinh))
                        return power(power(t,_ex2)-_ex1,_ex_1_2);

                // csch(asinh(x)) -> 1/x
                if (is_ex_the_function(x, asinh))
                        return power(t, _ex_1);
        }

        return csch(x).hold();
}

static ex csch_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx csch(x) -> -csch(x)*coth(x)
        return -mul(csch(x), coth(x));
}

static ex csch_series(const ex &x,
                      const relational &rel,
                      int order,
                      unsigned options)
{
        GINAC_ASSERT(is_a<symbol>(rel.lhs()));
        // method:
        // Taylor series where there is no pole falls back to tanh_deriv.
        // On a pole simply expand sinh(x)/cosh(x).
        const ex x_pt = x.subs(rel, subs_options::no_pattern);
        if (!(2*I*x_pt/Pi).info(info_flags::even))
                throw do_taylor();  // caught by function::series()
        // if we got here we have to care for a simple pole
        return (_ex1/sinh(x)).series(rel, order, options);
}

static ex csch_real_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return mul(cos(b), sinh(a)) / (power(mul(sin(b), cosh(a)), _ex2) + power(mul(cos(b), sinh(a)), _ex2));
}

static ex csch_imag_part(const ex & x)
{
        ex a = real_part(x);
        ex b = imag_part(x);
        return -mul(sin(b), cosh(a)) / (power(mul(sin(b), cosh(a)), _ex2) + power(mul(cos(b), sinh(a)), _ex2));
}

static ex csch_conjugate(const ex & x)
{
        return csch(x.conjugate());
}

REGISTER_FUNCTION(csch, eval_func(csch_eval).
                        evalf_func(csch_evalf).
                        derivative_func(csch_deriv).
                        series_func(csch_series).
                        real_part_func(csch_real_part).
                        imag_part_func(csch_imag_part).
                        conjugate_func(csch_conjugate).
                        latex_name("\\operatorname{csch}"));


//////////
// inverse hyperbolic cotangent (trigonometric function)
//////////

static ex acoth_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return atanh(ex_to<numeric>(x).inverse());

        return acoth(x).hold();
}

static ex acoth_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // acoth(1) -> oo
                if (x.is_equal(_ex1))
                        throw pole_error("acoth_eval(): simple pole",1);
                // acoth(-1) -> -oo
                if (x.is_equal(_ex_1))
                        throw pole_error("acoth_eval(): simple pole",1);
                //acoth(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return atanh(ex_to<numeric>(x).inverse());
                // acoth() is odd
                if (x.info(info_flags::negative))
                        return -acoth(-x);
        }

        if (is_exactly_a<GiNaC::function>(x)) {
                const ex &t = x.op(0);

                // acoth(coth(x)) -> x
                if (is_ex_the_function(x, coth))
                        return t;
        }


        return acoth(x).hold();
}

static ex acoth_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // acoth(x) -> (1/2)*(ln(1 + 1/x) - ln(1 - 1/x))
        // d/dx acoth(x) -> 1/(1-x^2)
        return power(_ex1-power(x, _ex2), _ex_1);
}

static ex acoth_conjugate(const ex & x)
{
        // conjugate(acoth(x))==acoth(conjugate(x)) unless on the branch cuts which
        // run along the real axis inside the interval [-1, +1].
        if (is_exactly_a<numeric>(x) &&
            (!x.imag_part().is_zero() || (x < *_num_1_p && x > *_num1_p))) {
                return acoth(x.conjugate());
        }
        return conjugate_function(acoth(x)).hold();
}

REGISTER_FUNCTION(acoth, eval_func(acoth_eval).
                         evalf_func(acoth_evalf).
                         derivative_func(acoth_deriv).
                         conjugate_func(acoth_conjugate).
                         set_name("arccoth"));

//////////
// inverse hyperbolic Cosecant (trigonometric function)
//////////

static ex acsch_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return asinh(ex_to<numeric>(x).inverse());

        return acsch(x).hold();
}

static ex acsch_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // acsch(0) -> oo
                if (x.is_zero())
                        throw pole_error("acsch_eval(): simple pole",1);
                //acsch(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return asinh(ex_to<numeric>(x).inverse());
                // acsch(-x) -> acsch(-x)
                if (x.info(info_flags::negative))
                        return -acsch(-x);
        }


        return acsch(x).hold();
}

static ex acsch_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // acsch(x) -> ln(1/x + sqrt(1/x^2 + 1))
        // d/dx acsch(x) ->  -1 / [x * sqrt(1 + x^2)];
        return (_ex_1/x)*power(_ex1+power(x, _ex2), _ex_1_2);
}

static ex acsch_conjugate(const ex & x)
{
        // conjugate(acsch(x))==acsch(conjugate(x)) unless on the branch cuts which
        // run along the imaginary axis inside the interval [-I, +I].
        if (x.info(info_flags::real))
                return acsch(x);
        if (is_exactly_a<numeric>(x)) {
                const numeric x_re = ex_to<numeric>(x.real_part());
                const numeric x_im = ex_to<numeric>(x.imag_part());
                if (!x_re.is_zero() ||
                    (x_im < *_num_1_p && x_im > *_num1_p))
                        return acsch(x.conjugate());
        }
        return conjugate_function(acsch(x)).hold();
}

REGISTER_FUNCTION(acsch, eval_func(acsch_eval).
                         evalf_func(acsch_evalf).
                         derivative_func(acsch_deriv).
                         conjugate_func(acsch_conjugate).
                         set_name("arccsch"));
//////////
// inverse hyperbolic Secant (trigonometric function)
//////////

static ex asech_evalf(const ex & x)
{
        if (is_exactly_a<numeric>(x))
                return acosh(ex_to<numeric>(x).inverse());

        return asech(x).hold();
}

static ex asech_eval(const ex & x)
{
        if (is_exactly_a<numeric>(x)) {
                // asech(0) -> oo
                if (x.is_zero())
                        throw pole_error("asech_eval(): simple pole",1);
                // asech(1) -> 0
                if (x.is_equal(_ex1))
                        return _ex0;
                //asech(-1) -> I*Pi
                if (x.is_equal(_ex_1))
                        return Pi*I;
                //asech(float) -> float
                if (x.info(info_flags::numeric) && !x.info(info_flags::crational))
                        return acosh(ex_to<numeric>(x).inverse());
                // asech(-x) -> Pi*I-asech(-x)
                if (x.info(info_flags::negative))
                        return Pi*I-asech(-x);
        }


        return asech(x).hold();
}

static ex asech_deriv(const ex & x, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // asech(x) -> ln(1/x + sqrt(1/x^2 - 1))
        // d/dx asech(x) ->  -1 / [x * sqrt(1 - x^2)];
        return (_ex_1/x)*power(_ex1-power(x, _ex2), _ex_1_2);
}

static ex asech_conjugate(const ex & x)
{
        // conjugate(asech(x))==asech(conjugate(x)) unless on the branch cuts which
        // run along the real axis from 0 to -oo and 1 to oo.
        if (is_exactly_a<numeric>(x) &&
            (!x.imag_part().is_zero() || (x < *_num1_p && x > *_num0_p))) {
                return asech(x.conjugate());
        }
        return conjugate_function(asech(x)).hold();
}

REGISTER_FUNCTION(asech, eval_func(asech_eval).
                         evalf_func(asech_evalf).
                         derivative_func(asech_deriv).
                         conjugate_func(asech_conjugate).
                         set_name("arcsech"));


///////////
/// JacobiSN function
///////////

static ex JacobiSN_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiSN(x,m) ->  JacobiCN(x, m)*JacobiDN(x, m);
        return JacobiCN(x, m)*JacobiDN(x, m);
}

REGISTER_FUNCTION(JacobiSN, derivative_func(JacobiSN_deriv))


///////////
/// JacobiCN functions
///////////

static ex JacobiCN_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiCN(x,m) ->  -JacobiDN(x, m)*JacobiSN(x, m);
        return -JacobiDN(x, m)*JacobiSN(x, m);
}

REGISTER_FUNCTION(JacobiCN, derivative_func(JacobiCN_deriv))


///////////
/// JacobiDN functions
///////////

static ex JacobiDN_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiDN(x,m) -> -m^2*JacobiCN(x, m)*JacobiSN(x, m);
        return -m*m*JacobiCN(x, m)*JacobiSN(x, m);
}

REGISTER_FUNCTION(JacobiDN, derivative_func(JacobiDN_deriv))


///////////
/// JacobiNS functions
///////////

static ex JacobiNS_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiNS(x,m) -> -JacobiCN(x, m)*JacobiDN(x, m)/JacobiSN(x, m)^2;
        return -JacobiCN(x, m)*JacobiDN(x, m)/pow(JacobiSN(x, m),2);
}

REGISTER_FUNCTION(JacobiNS, derivative_func(JacobiNS_deriv))


///////////
/// JacobiNC functions
///////////

static ex JacobiNC_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiNC(x,m) -> JacobiDN(x, m)*JacobiSN(x, m)/JacobiCN(x, m)^2;
        return JacobiDN(x, m)*JacobiSN(x, m)/pow(JacobiCN(x, m),2);
}

REGISTER_FUNCTION(JacobiNC, derivative_func(JacobiNC_deriv))


///////////
/// JacobiND functions
///////////

static ex JacobiND_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiND(x,m) -> m^2*JacobiCN(x, m)*JacobiSN(x, m)/JacobiDN(x, m)^2;
        return m*m*JacobiCN(x, m)*JacobiSN(x, m)/pow(JacobiDN(x, m),2);
}

REGISTER_FUNCTION(JacobiND, derivative_func(JacobiND_deriv))


///////////
/// JacobiSC functions
///////////

static ex JacobiSC_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiSC(x,m) -> JacobiDN(x, m)/JacobiCN(x, m)^2;
        return JacobiDN(x, m)/pow(JacobiCN(x, m),2);
}

REGISTER_FUNCTION(JacobiSC, derivative_func(JacobiSC_deriv))


///////////
/// JacobiSD functions
///////////

static ex JacobiSD_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiSD(x,m) -> JacobiCN(x, m)/JacobiDN(x, m)^2;
        return JacobiCN(x, m)/pow(JacobiDN(x, m),2);
}

REGISTER_FUNCTION(JacobiSD, derivative_func(JacobiSD_deriv))


///////////
/// JacobiNC functions
///////////

static ex JacobiCS_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiCS(x,m) -> -JacobiDN(x, m)/JacobiSN(x, m)^2;
        return -JacobiDN(x, m)/pow(JacobiSN(x, m),2);
}

REGISTER_FUNCTION(JacobiCS, derivative_func(JacobiCS_deriv))


///////////
/// JacobiDS functions
///////////

static ex JacobiDS_deriv(const ex & x,const ex & m, unsigned deriv_param)
{
        GINAC_ASSERT(deriv_param==0);

        // d/dx JacobiDS(x,m) -> -JacobiCS(x, m)/JacobiSN(x, m);
        return -JacobiCS(x, m)/JacobiSN(x, m);
}

REGISTER_FUNCTION(JacobiDS, derivative_func(JacobiDS_deriv))




//////////
// Derivative of Diff function
//////////

static ex Diff_deriv(const ex& expr_, const ex& var_1, const ex& order_1, unsigned diff_param)
{
    return Diff(expr_, var_1, order_1 + _ex1);

}

REGISTER_FUNCTION(Diff,derivative_func(Diff_deriv));



//////////
// Integration of Integrate function
//////////

static ex Integrate_deriv(const ex& expr_, const ex& var_1, const ex& order_1, unsigned diff_param)
{
    if (order_1 > _ex0)
        return Integrate(expr_, var_1, order_1 + _ex_1);
    else if (order_1 == _ex0)
        return Diff(expr_, var_1, _ex1);
    else
        return _ex0;
}

REGISTER_FUNCTION(Integrate, derivative_func(Integrate_deriv));


///////
// generators of all symbols and expressions.
// declare here to include all above functions.
//////
parser reader;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
