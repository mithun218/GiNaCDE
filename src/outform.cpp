
/** @file outform.cpp
 *
 *  Implementation of GiNaCDE's output formt. Output results can be saved
    in Maple or Mathematica format. */



#include <string.h>
#include <string>
#include <cctype>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include "outform.h"
#include "desolve.h"
#include "utility.h"

using namespace std;

int output = maple;
string dpndtWtIndpndt = "u";

/// output structure ///
string outstr(const char* _sym, int symno)
{
    char temstr[symno];

    strcpy(temstr, _sym);
    for(int i = 0 ; i< symno - 1; i++)
        strcat(temstr, _sym);

    string str = temstr;

    return str;
}

/// replacing substring in string ///
string replacestring(std::string subject, const string& search,
const string& replace)
{
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos)
    {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

/// Bracket matcing ///
bool bktmch(const string _instr)
{
    string obkt = "({[", cbkt = ")}]";
    map<char, char> cobkt;
    vector<char> store;

    cobkt[')'] = '(';
    cobkt['}'] = '{';
    cobkt[']'] = '{';

    for(string::const_iterator it = _instr.begin(); it != _instr.end(); it++)
    {
        if(obkt.find(*it) != string::npos)
            store.push_back(*it);
        if(cbkt.find(*it) != string::npos)
        {
            if(store.empty())
                return false;

            char tem = store.back();
            if(tem == cobkt[*it])
                store.pop_back();
        }
    }

    if(store.empty())
        return true;
    else
        return false;
}

/// ginac output format to mathematica output format ///
string gmathematica(string _instr)
{
    size_t instrlen = _instr.length();

    string obkt = "({[", cbkt = ")}]";
    map<char, char> cobkt;
    vector<char> store;

    cobkt[')'] = '(';
    cobkt['}'] = '{';
    cobkt[']'] = '[';

    bool cond1 = true;

    size_t strindx = 0;
    size_t prevstrindx = 0;

    while(cond1)
    {
        if(strindx < instrlen && isalpha(_instr[strindx]))
        {
             //passing  alphabets//
             if(strindx + 1 < instrlen && isalpha(_instr[strindx + 1]))
                strindx = strindx + 1;

             else if (strindx + 1 < instrlen && _instr[strindx + 1] == '(')
             {
                _instr[prevstrindx] = toupper(_instr[prevstrindx]);

                // searching next closing bracket ')'//
                bool cond2 = true;
                size_t cbktpos = ++strindx;
                while(cond2)
                {
                    if(obkt.find(_instr[cbktpos]) != string::npos)
                        store.push_back(_instr[cbktpos]);
                    if(cbkt.find(_instr[cbktpos]) != string::npos)
                    {
                        char tem = store.back();
                        if(tem == cobkt[_instr[cbktpos]])
                            store.pop_back();
                    }

                    if(store.empty())
                        cond2 = false;
                    else
                        cbktpos = cbktpos + 1;
                }
                _instr[strindx] = '[';
                _instr[cbktpos] = ']';

             }

             else
                strindx = strindx + 1;

        }

        else if (strindx >= instrlen)
            cond1 = false;
        else
        {
            strindx = strindx + 1;
            prevstrindx = strindx;
        }

    }

    return _instr;
}

/** converting output form of input equation into maple or mathematica form **/
string diffformchange(const ex& diffeq, const lst& dpndt_vars, const exset& indpndt_vars)
{

    stringstream tem1, tem2, solutionsdiffst, dpndt_varsStr;
    string solutionsdiff;

    solutionsdiffst << diffeq;
    solutionsdiff = solutionsdiffst.str();

    if(output==ginac)
        return solutionsdiffst.str();

    for(size_t dpndti=0;dpndti<nops(dpndt_vars);dpndti++)
    {
        dpndt_varsStr << dpndt_vars.op(dpndti);
        const string dpndt_varsStr2 = dpndt_varsStr.str();

        if(output == maple)
        {
            tem2 << dpndt_varsStr2 << "(";
            for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
            {
                tem2 << *it;
                if(next(it) != indpndt_vars.end())
                    tem2 << ",";
            }
            tem2 << ")";
        }
        else if( output == mathematica )
        {
            tem2 << dpndt_varsStr2 << "[";

            for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
            {
                tem2 << *it;
                if(next(it) != indpndt_vars.end())
                    tem2 << ",";
            }
            tem2 << "]";
        }

        if(dpndtWtIndpndt == "u")
            dpndtWtIndpndt = tem2.str();

        vector< string > tem1Clt,tem2Clt;

        tem1Clt.push_back( "+" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "+"+tem2.str()+"+" );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "-"+tem2.str()+"+" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "*"+tem2.str()+"+" );
        tem1Clt.push_back( "(" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "("+tem2.str()+"+" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "^"+tem2.str()+"+" );
        tem1Clt.push_back( "(" + dpndt_varsStr2 + "+" );
        tem2Clt.push_back( "("+tem2.str()+"+" );
        tem1Clt.push_back("[" + dpndt_varsStr2 + "+");
        tem2Clt.push_back( "["+tem2.str()+"+" );


        tem1Clt.push_back( "+" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "+"+tem2.str()+"-" );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "-"+tem2.str()+"-" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "*"+tem2.str()+"-" );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "/"+tem2.str()+"-" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "^"+tem2.str()+"-" );
        tem1Clt.push_back("(" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "("+tem2.str()+"-" );
        tem1Clt.push_back( "[" + dpndt_varsStr2 + "-" );
        tem2Clt.push_back( "["+tem2.str()+"-" );


        tem1Clt.push_back( "+" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "+"+tem2.str()+"*" );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "-"+tem2.str()+"*" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "*"+tem2.str()+"*" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "*"+tem2.str()+"*" );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "/"+tem2.str()+"*" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "^"+tem2.str()+"*" );
        tem1Clt.push_back( "(" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "("+tem2.str()+"*" );
        tem1Clt.push_back( "[" + dpndt_varsStr2 + "*" );
        tem2Clt.push_back( "["+tem2.str()+"*" );


        tem1Clt.push_back( "+" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "+"+tem2.str()+"/" );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "/"+tem2.str()+"/" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "*"+tem2.str()+"/" );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "/"+tem2.str()+"/" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "^"+tem2.str()+"/" );
        tem1Clt.push_back( "(" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "("+tem2.str()+"/" );
        tem1Clt.push_back( "[" + dpndt_varsStr2 + "/" );
        tem2Clt.push_back( "["+tem2.str()+"/" );


        tem1Clt.push_back( "+" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "+"+tem2.str()+"^" );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "-"+tem2.str()+"^" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "^"+tem2.str()+"^" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "*"+tem2.str()+"^" );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "/"+tem2.str()+"^" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "^"+tem2.str()+"^" );
        tem1Clt.push_back( "(" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "("+tem2.str()+"^" );
        tem1Clt.push_back( "[" + dpndt_varsStr2 + "^" );
        tem2Clt.push_back( "["+tem2.str()+"^" );


        tem1Clt.push_back( "(" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "("+tem2.str()+"," );

        tem1Clt.push_back( "," + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( ","+tem2.str()+")" );
        tem1Clt.push_back( "+" + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( "+"+tem2.str()+")" );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( "-"+tem2.str()+")" );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( "*"+tem2.str()+")" );
        tem1Clt.push_back( "^" + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( "^"+tem2.str()+")" );




        tem1Clt.push_back( "[" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "["+tem2.str()+"," );

        tem1Clt.push_back( "," + dpndt_varsStr2 + "]" );
        tem2Clt.push_back( ","+tem2.str()+"]" );


        tem1Clt.push_back( "[" + dpndt_varsStr2 + "]" );
        tem2Clt.push_back( "["+tem2.str()+"]" );

        tem1Clt.push_back( "(" + dpndt_varsStr2 + ")" );
        tem2Clt.push_back( "("+tem2.str()+")" );


        tem1Clt.push_back( "," + dpndt_varsStr2 + "," );
        tem2Clt.push_back( ","+tem2.str()+"," );
        tem1Clt.push_back( "*" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "*"+tem2.str()+"," );
        tem1Clt.push_back( "+" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "+"+tem2.str()+"," );
        tem1Clt.push_back( "-" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "-"+tem2.str()+"," );
        tem1Clt.push_back( "/" + dpndt_varsStr2 + "," );
        tem2Clt.push_back( "/"+tem2.str()+"," );

        const size_t strleng = dpndt_varsStr2.length();
        if(solutionsdiff.length()> strleng)
        {
            // Taking first and last two variables
            const string temstr1 = solutionsdiff.substr(0,strleng+1);
            const string temstr2 = solutionsdiff.substr(solutionsdiff.length()-strleng-1);

            // changing in first dependent variable
            if(temstr1.substr(0,strleng) == dpndt_varsStr2 && (temstr1[strleng]=='+' || temstr1[strleng]=='-' || temstr1[strleng]=='*' || temstr1[strleng]=='/' || temstr1[strleng]=='^'))
            {
                solutionsdiff.replace(0, strleng, tem2.str());
            }
            // changing in last dependent variable
            else if(temstr2.substr(1) == dpndt_varsStr2 && (temstr2[0]=='+' || temstr2[0]=='-' || temstr2[0]=='*' || temstr2[0]=='/' || temstr2[0]=='^'))
            {
                solutionsdiff.replace(solutionsdiff.length()-strleng, strleng, tem2.str());
            }
        }

        if( output == maple )
        {
            for( unsigned i = 0; i < tem1Clt.size(); i++ )
            {
                solutionsdiff = replacestring(solutionsdiff, tem1Clt[i], tem2Clt[i]);
            }

            solutionsdiff = (replacestring(solutionsdiff, "Diff", "diff"));
        }
        else if( output == mathematica )
        {

            for( unsigned i = 0; i < tem1Clt.size(); i++ )
            {
                solutionsdiff = replacestring(solutionsdiff, tem1Clt[i], tem2Clt[i]);
            }

            solutionsdiff = replacestring(solutionsdiff, "Diff", "D");
        }


        tem1.str("");
        tem2.str("");
        dpndt_varsStr.str("");
        tem1Clt.clear();
        tem2Clt.clear();
    }


    if(output==maple)
    {
        tem1.str("");
        tem2.str("");

        int difford = order(diffeq);
        for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
        {
            int varcount = 1;
            do
            {
                tem1 << *it << "," << varcount;
                tem2 << *it << "$" << varcount;
                varcount++;
                solutionsdiff = (replacestring(solutionsdiff, tem1.str(), tem2.str()));
                tem1.str("");
                tem2.str("");
            }while(varcount <= difford);
        }

    }
    else if(output==mathematica)
    {
        tem1.str("");
        tem2.str("");
        int difford = order(diffeq);
        for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
        {
            int varcount = 1;
            do
            {
                tem1 << *it << "," << varcount;
                tem2 << "{" << *it << "," << varcount << "}";
                varcount++;
                solutionsdiff = (replacestring(solutionsdiff, tem1.str(), tem2.str()));
                tem1.str("");
                tem2.str("");
            }while(varcount <= difford);
        }

    }

    return solutionsdiff;
}
/////////////////////////////////////////////////////////////


string writetofile(stringstream& stringbuf, const ex& dpndt_var)
{    
    ofstream outfile;
    stringstream dpndt_varStr, solutionStr;

    dpndt_varStr << dpndt_var<< " =";

    if(output == mathematica)
    {
        solutionStr << replacestring(gmathematica(replacestring(replacestring(replacestring(replacestring(replacestring(replacestring(replacestring(stringbuf.str(),
                   "g_", "gun"), "h_", "hun"), "Y_", "Yun"), "X_", "Xun"), "C_", "Const"),"==", "="),"_","")),dpndt_varStr.str(),dpndtWtIndpndt+" =");

        outfile.open(filename);
        outfile<<solutionStr.str();
        outfile.close();
    }
    else if(output == maple)
    {
        solutionStr << replacestring(replacestring(stringbuf.str(), "==", "="), dpndt_varStr.str(), dpndtWtIndpndt+" =");

        outfile.open(filename);
        outfile<<solutionStr.str();
        outfile.close();
    }
    else
    {
        outfile.open(filename);
        outfile << stringbuf.str();
        outfile.close();

        return stringbuf.str();
    }

    dpndtWtIndpndt = "u";

    return solutionStr.str();

}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////























