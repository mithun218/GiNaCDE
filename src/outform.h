
/** @file outform.h
 *
 *  Interface to the GiNaCDE's output formt implemented in outform.cpp. */



#ifndef OUTFORM_H_INCLUDED
#define OUTFORM_H_INCLUDED

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

#define maple 8
#define mathematica 9
#define ginac 10

extern int output;
extern string dpndtWtIndpndt; // It is assinged to dependend variable with independen variable(s), such as u(t,x)

string outstr(const char* _sym, int symno);

string replacestring(std::string subject, const string& search,
const string& replace);

bool bktmch(const string _instr);

string gmathematica(string _instr);

string diffformchange(const ex& diffeq, const lst& dpndt_var, const exset& indpndt_var);

void writetofile(stringstream&, const ex& dpndt_var);

#endif // OUTFORM_H_INCLUDED
