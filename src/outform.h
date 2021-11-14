
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

string outstr(const char* _sym, int symno);

string replacestring(std::string subject, const string& search,
const string& replace);

bool bktmch(const string _instr);

string gmathematica(string _instr);

string diffformchange(const ex& diffeq, const lst& dpndt_var, const exset& indpndt_var);

void writetofile(stringstream&);

#endif // OUTFORM_H_INCLUDED
