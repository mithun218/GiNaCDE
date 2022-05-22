
/** @file inifcns.h
 *
 *  Interface to some usefull funtions implemented in inifcns.cpp. */



#ifndef INIFCNS_H_INCLUDED
#define INIFCNS_H_INCLUDED

#define Real 1
#define posreal 2
#define noassume 3

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;


/** Inverse cotangent (arc cotangent). */
DECLARE_FUNCTION_1P(acot)

/** Inverse secant (arc secant). */
DECLARE_FUNCTION_1P(asec)

/** Inverse cosecant (arc cosecant). */
DECLARE_FUNCTION_1P(acsc)

/** Cotangent. */
DECLARE_FUNCTION_1P(cot)

/** Secant. */
DECLARE_FUNCTION_1P(sec)

/** Cosecant. */
DECLARE_FUNCTION_1P(csc)


/** Hyperbolic Cotangent. */
DECLARE_FUNCTION_1P(coth)

/** Hyperbolic Secant. */
DECLARE_FUNCTION_1P(sech)

/** Hyperbolic Cosecant. */
DECLARE_FUNCTION_1P(csch)

/** Inverse hyperbolic Cotangent (area hyperbolic cotangent). */
DECLARE_FUNCTION_1P(acoth)

/** Inverse hyperbolic Cosecant (area hyperbolic cosecant). */
DECLARE_FUNCTION_1P(acsch)

/** Inverse hyperbolic Secant (area hyperbolic secant). */
DECLARE_FUNCTION_1P(asech)

/** Dummy Jacobi Elliptic Functions **/
DECLARE_FUNCTION_2P(JacobiSN)
DECLARE_FUNCTION_2P(JacobiCN)
DECLARE_FUNCTION_2P(JacobiDN)
DECLARE_FUNCTION_2P(JacobiNS)
DECLARE_FUNCTION_2P(JacobiNC)
DECLARE_FUNCTION_2P(JacobiND)
DECLARE_FUNCTION_2P(JacobiSC)
DECLARE_FUNCTION_2P(JacobiSD)
DECLARE_FUNCTION_2P(JacobiCS)
DECLARE_FUNCTION_2P(JacobiDS)


/** Noun form diff for one variable **/
DECLARE_FUNCTION_3P(Diff)

/** Noun forn integrate **/
DECLARE_FUNCTION_3P(Integrate)

/** generator of expressions **/
extern parser reader;
//extern map<string, ex> directory;
//extern ex reader(const string& _expr, const int& assume=posreal);

#endif // INIFCNS_H_INCLUDED
