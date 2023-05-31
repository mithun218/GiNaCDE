
/** @file desolve.cpp
 *
 *  Implementation of the function desolve*/



#include <ginac/ginac.h>
#include <sstream>
#include <chrono>
#include "utility.h"
#include "derivative.h"
#include "integrate.h"
#include "solve.h"
#include "F_expns_methd.h"
#include "desolve.h"
#include "outform.h"
#include "simplify.h"
#include "fim.h"
#include "inifcns.h"


using namespace std;
using namespace GiNaC;

ex chi, xi;

ex NValue = 0;

lst twcPhase= {}, degAcoeff={}, paraInDiffSolve = {};





#ifdef GiNaCDE_gui
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

vector<GtkWidget*> entry, entrylbl;
GtkWidget *phaseCombo;
vector<ex> entryText;

bool cancel = false;


static int phasepart_okbutton_clicked(GtkWidget* okbutton,gpointer dialog)
{
    //GdkColor gtk_color;
    //gdk_color_parse ("#FF0000", &gtk_color);

    bool isValid = true;
    entryText.clear();
    for(unsigned i=0; i<entry.size(); i++)
    {
        try
        {
            entryText.push_back(reader(gtk_entry_get_text(GTK_ENTRY(entry[i]))));
        }
        catch(GiNaC::parse_error)
        {
            isValid = false;

            const gchar *strLbl = g_strdup_printf ("<span color=\"red\">"
               "%s"
             "</span>",
             gtk_label_get_text(GTK_LABEL(entrylbl[i])));
            gtk_label_set_markup(GTK_LABEL(entrylbl[i]), strLbl);

            //gtk_widget_modify_fg(entrylbl[i], GTK_STATE_NORMAL, &gtk_color);
            
            g_free((gpointer) strLbl);
            
        }

    }

    if(!isValid)
        return 0;

    gtk_widget_destroy(GTK_WIDGET(dialog));
    return 0;
}

static gboolean key_press_cb(GtkWidget *w, GdkEvent *ev, gpointer data)
{
    GdkEventKey *key = (GdkEventKey*)ev;
    if(key->keyval == GDK_KEY_Escape) // preventing closing of dialog for escape key
    {
        return TRUE;
    }

    const gchar *strLbl = g_strdup_printf ("<span color=\"red\">"
       "%s"
     "</span>",
     gtk_label_get_text(GTK_LABEL(data)));
    gtk_label_set_markup(GTK_LABEL(data), strLbl);

    //gtk_widget_modify_fg(GTK_WIDGET(data), GTK_STATE_NORMAL, NULL);
    g_free((gpointer) strLbl);
    return FALSE;
}

static int twc_okbutton_clicked(GtkWidget* okbutton,gpointer dialog)
{
    bool isValid = true;
    //GdkColor gtk_color;
    //gdk_color_parse ("#FF0000", &gtk_color);

    entryText.clear();
    for(unsigned i=0; i<entry.size(); i++)
    {
        try
        {
            entryText.push_back(reader(gtk_entry_get_text(GTK_ENTRY(entry[i]))));
        }
        catch(GiNaC::parse_error)
        {
            isValid = false;

            const gchar *strLbl = g_strdup_printf ("<span color=\"red\">"
               "%s"
             "</span>",
            gtk_label_get_text(GTK_LABEL(entrylbl[i])));
            gtk_label_set_markup(GTK_LABEL(entrylbl[i]), strLbl);

            //gtk_widget_modify_fg(entrylbl[i], GTK_STATE_NORMAL, &gtk_color);
            g_free((gpointer) strLbl);
        }

    }

    if(!isValid)
        return 0;

   gtk_widget_destroy(GTK_WIDGET(dialog));
   return 0;
}


static void cancelbutton_clicked(GtkWidget* cancelbutton,gpointer dialog)
{
    cancel = true;
    gtk_widget_destroy(GTK_WIDGET(dialog));
    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
}

#endif // GiNaCDE_gui

int odeType_check(const ex& _ode, const ex& _dpndtVar)
{
    exmap repls2;
    const ex C1_=reader("C1_");
    int odetype = general;

    (Collect(_ode, _dpndtVar)).match( wild(2)*pow(_dpndtVar, wild(3)) + wild(1)*_dpndtVar, repls2);
    if(!repls2.empty()&& !repls2[wild(1)].has(_dpndtVar)&& !repls2[wild(2)].has(_dpndtVar)&& !repls2[wild(3)].has(_dpndtVar) )
        return bernouli;

    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( wild(2)*pow(_dpndtVar, wild(3)), repls2);
    if(!repls2.empty()&& !repls2[wild(2)].has(_dpndtVar)&& !repls2[wild(3)].has(_dpndtVar))
    {
        if((repls2[wild(3)]).info(info_flags::numeric) && repls2[wild(3)].info(info_flags::posint))
            return bernouli;
        else
            return fracbernouli; // handling y' = A*y^(n/2)
    }

    // y' = sqrt(A*y)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match(sqrt( wild(2)*_dpndtVar), repls2);
    if(!repls2.empty()&& !repls2[wild(2)].has(_dpndtVar))
    {
            return fracbernouli;
    }


    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( wild(1)*_dpndtVar, repls2);
    if(!repls2.empty()&& !repls2[wild(1)].has(_dpndtVar))
        return bernouli;

//////////////////////////////////////////////////////////////////////
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( wild(2)*pow(_dpndtVar, _ex2) + wild(1)*_dpndtVar + wild(0), repls2);
    if(!repls2.empty() && !repls2[wild(2)].has(_dpndtVar) && !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(0)].has(_dpndtVar))
        return riccati;

    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( wild(2)*pow(_dpndtVar, _ex2) + wild(0), repls2);
    if(!repls2.empty() && !repls2[wild(2)].has(_dpndtVar) && !repls2[wild(0)].has(_dpndtVar))
        return riccati;

    repls2.clear();
    (Collect(_ode, _dpndtVar)).match(  wild(1)*_dpndtVar + wild(0), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar) && !repls2[wild(1)].has(_dpndtVar))
    {
        return riccati;
    }

    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( wild(0), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar))
    {
        return riccati;
    }


///////////////////////////////////////////////////////////////////
// sqrt(a2*y(t)^2+a4*y(t)^4)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(2)*pow(_dpndtVar, 4)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype1;
    }

///////////////////////////////////////////////////////////////////
// sqrt(a0+a2*y(t)^2)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(2)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype2;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a2*y(t)^2+a3*y(t)^3)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(2)*pow(_dpndtVar, 3)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype3;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a2*y(t)^2+a4*y(t)^4+a6*y(t)^6)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(1)*pow(_dpndtVar,4) + wild(2)*pow(_dpndtVar,6)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype4;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a2*y(t)^2+a6*y(t)^6)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(2)*pow(_dpndtVar,6)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype4;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a0+a2*y(t)^2+a4*y(t)^4+a6*y(t)^6)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(3) + wild(0)*pow(_dpndtVar, 2) + wild(1)*pow(_dpndtVar,4) + wild(2)*pow(_dpndtVar,6)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar)
       && !repls2[wild(3)].has(_dpndtVar))
    {
        return odetype5;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a2*y(t)^2+a3*y(t)^3+a4*y(t)^4)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*pow(_dpndtVar, 2) + wild(1)*pow(_dpndtVar,3) + wild(2)*pow(_dpndtVar,4)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar))
    {
        return odetype6;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a0+a1*y(t))
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)+ wild(1)*_dpndtVar), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar))
    {
        return odetype7;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a0+a2*y(t)^2+a4*y(t)^4)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0) + wild(1)*pow(_dpndtVar,2) + wild(2)*pow(_dpndtVar,4)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar))
    {
        return jacobiElip024;
    }

//////////////////////////////////////////////////////////////////
// sqrt(a1*y(t)+a2*y(t)^2+a3*y(t)^3)
    repls2.clear();
    (Collect(_ode, _dpndtVar)).match( sqrt(wild(0)*_dpndtVar + wild(1)*pow(_dpndtVar,2) + wild(2)*pow(_dpndtVar,3)), repls2);
    if(!repls2.empty() && !repls2[wild(0)].has(_dpndtVar)&& !repls2[wild(1)].has(_dpndtVar) && !repls2[wild(2)].has(_dpndtVar))
    {
        return jacobiElip123;
    }

    return odetype;
}

/** returning solution for Bernouli, Riccati diff equations **/
lst firstOrderDiff_solu(const ex& _indpndt_var, const int& _odetype)
{
    const ex C_=reader("C_");

    if(_odetype == riccati)
    {
        if(degAcoeff[0]==2&&( !(degAcoeff[3] == _ex0) && (( !(degAcoeff[2] == _ex0)) || ( !(degAcoeff[1] == _ex0))) ))
        {
            const ex theta = sqrt(simplify(degAcoeff[2]*degAcoeff[2]-_ex2*_ex2*degAcoeff[1]*degAcoeff[3])),
                     A1byA2 = simplify(-_ex1_2*degAcoeff[2]/degAcoeff[3]),
                     A1byA0 = simplify(-_ex1_2*degAcoeff[2]/degAcoeff[1]);

            // Checking theta is real //
            if (is_a<numeric>(theta) && !theta.info(info_flags::real))
            {
                return {};
            }

            return {A1byA2 - theta * tanh(_ex1_2 * theta * _indpndt_var + C_) / (2 * degAcoeff[3]),
                    A1byA2 - theta * coth(_ex1_2 * theta * _indpndt_var + C_) / (2 * degAcoeff[3]),
                    1
                        / (A1byA0
                           + theta * tanh((_ex1_2) *theta * _indpndt_var + C_) / (2 * degAcoeff[1])),
                    1
                        / (A1byA0
                           + theta * coth((_ex1_2) *theta * _indpndt_var + C_) / (2 * degAcoeff[1])),
                    A1byA2 - theta * tanh((_ex1_2) *theta * _indpndt_var + C_) / (2 * degAcoeff[3])
                        + pow(sech(_ex1_2 * theta * _indpndt_var + C_),_ex2)
                              / (C_
                                 - 2 * degAcoeff[3] * tanh((_ex1_2) *theta * _indpndt_var + C_)
                                       / theta)};
        }
        else if(degAcoeff[0]==1 && !(degAcoeff[2] == _ex0))
            return {simplify(-degAcoeff[1]/degAcoeff[2])+exp(degAcoeff[2]*_indpndt_var)*C_};
        else
            return {degAcoeff[1]*_indpndt_var +C_};
    }
    else if(_odetype == bernouli)
    {
        const ex odedeg = degAcoeff[0];
        if( odedeg > _ex1 )
        {
            if( degAcoeff[2] != _ex0 )
            {
                return {pow((degAcoeff[2]*(cosh(degAcoeff[2]*_indpndt_var*(odedeg-1)+C_*degAcoeff[2])+sinh(degAcoeff[2]*_indpndt_var*(odedeg-1)+C_*degAcoeff[2]))/(1-degAcoeff[odedeg+1]*cosh(degAcoeff[2]*_indpndt_var*(odedeg-1)+C_*degAcoeff[2])-degAcoeff[odedeg+1]*sinh(degAcoeff[2]*_indpndt_var*(odedeg-1)+C_*degAcoeff[2]))),(1/(odedeg-1))),
                       pow((-degAcoeff[odedeg+1]/degAcoeff[2]+C_*exp(degAcoeff[2]*(1-odedeg)*_indpndt_var)),(1/(1-odedeg))),
                       pow((-degAcoeff[2]/(2*degAcoeff[odedeg+1])-degAcoeff[2]*tanh((_ex1_2)*(odedeg-1)*degAcoeff[2]*_indpndt_var+C_)/(2*degAcoeff[odedeg+1])),(1/(odedeg-1))),
                       pow((-degAcoeff[2]/(2*degAcoeff[odedeg+1])-degAcoeff[2]*coth((_ex1_2)*(odedeg-1)*degAcoeff[2]*_indpndt_var+C_)/(2*degAcoeff[odedeg+1])),(1/(odedeg-1)))
                       };
            }
            else
            {
                return {_ex1/pow((_indpndt_var*degAcoeff[odedeg+1]*(_ex1-odedeg)+C_),(_ex1/(odedeg-1)))};
            }
        }
        else
          return {C_*exp( (degAcoeff[2]*_indpndt_var) )};

    }
    else if(_odetype == fracbernouli)
    {
        return {_ex1/pow((_indpndt_var*degAcoeff[1]*(_ex1-degAcoeff[0])+C_),(_ex1/(degAcoeff[0]-1)))};
    }
    else if(_odetype == odetype1)
    {
        const ex A2A4 = simplify(degAcoeff[3]*degAcoeff[5]);
        return {4*degAcoeff[3]*exp((C_+_indpndt_var)*sqrt(degAcoeff[3]))/(-4*A2A4*exp(2*_indpndt_var*sqrt(degAcoeff[3]))+exp(2*C_*sqrt(degAcoeff[3]))),-4*degAcoeff[3]*exp((C_+_indpndt_var)*sqrt(degAcoeff[3]))/(4*A2A4*exp(2*C_*sqrt(degAcoeff[3]))-exp(2*_indpndt_var*sqrt(degAcoeff[3]))),
        sqrt(-A2A4)*sech(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[5],-sqrt(-A2A4)*sech(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[5],sqrt(A2A4)*csch(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[5],-sqrt(A2A4)*csch(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[5]};
    }
    else if(_odetype == odetype2)
    {
       const ex A0A2 = simplify(degAcoeff[3]*degAcoeff[1]);
       return {-(_ex1_2)*(degAcoeff[1]*exp(2*C_*sqrt(degAcoeff[3]))-exp(2*_indpndt_var*sqrt(degAcoeff[3])))*exp(-(C_+_indpndt_var)*sqrt(degAcoeff[3]))/sqrt(degAcoeff[3]),(_ex1_2)*(-degAcoeff[1]*exp(2*_indpndt_var*sqrt(degAcoeff[3]))+exp(2*C_*sqrt(degAcoeff[3])))*exp(-(C_+_indpndt_var)*sqrt(degAcoeff[3]))/sqrt(degAcoeff[3]),
       sqrt(-A0A2)*cosh(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[3],-sqrt(-A0A2)*cosh(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[3],sqrt(A0A2)*sinh(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[3],-sqrt(A0A2)*sinh(_indpndt_var*sqrt(degAcoeff[3])+C_)/degAcoeff[3]};
    }
    else if(_odetype == odetype3)
    {
        const ex A2byA3 = degAcoeff[3]/degAcoeff[4];
        return {-A2byA3*pow(sech((_ex1_2)*_indpndt_var*sqrt(degAcoeff[3])),2),A2byA3*pow(csc((_ex1_2)*_indpndt_var*sqrt(degAcoeff[3])),2)};
    }
    else if(_odetype == odetype4)
    {
        const exset symClt = symbols(degAcoeff[3]+degAcoeff[5]+degAcoeff[7]);

        if(!symClt.empty()) // minimum one symbol is required for two additional solutions with conditions, a6=a4^2/(4*a2)
          {
              const ex theta = sqrt(simplify(degAcoeff[5]*degAcoeff[5]-_ex2*_ex2*degAcoeff[3]*degAcoeff[7]));
              return {pow((2*degAcoeff[3]*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta-(theta+degAcoeff[5])*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
                     -pow((2*degAcoeff[3]*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta-(theta+degAcoeff[5])*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
                      pow((2*degAcoeff[3]*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta+(theta-degAcoeff[5])*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
                     -pow((2*degAcoeff[3]*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta+(theta-degAcoeff[5])*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
                      lst{pow((-degAcoeff[3]*(1+tanh(sqrt(degAcoeff[3])*_indpndt_var+C_))/degAcoeff[5]),(_ex1_2)),degAcoeff[7]==pow(degAcoeff[5],2)/(4*degAcoeff[3])},
                      lst{pow((-degAcoeff[3]*(1+coth(sqrt(degAcoeff[3])*_indpndt_var+C_))/degAcoeff[5]),(_ex1_2)),degAcoeff[7]==pow(degAcoeff[5],2)/(4*degAcoeff[3])}
                      };
          }

        const ex theta = sqrt(simplify(degAcoeff[5]*degAcoeff[5]-_ex2*_ex2*degAcoeff[3]*degAcoeff[7]));
        return {pow((2*degAcoeff[3]*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta-(theta+degAcoeff[5])*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
               -pow((2*degAcoeff[3]*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta-(theta+degAcoeff[5])*pow(sech(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
                pow((2*degAcoeff[3]*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta+(theta-degAcoeff[5])*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2),
               -pow((2*degAcoeff[3]*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2)/(2*theta+(theta-degAcoeff[5])*pow(csch(sqrt(degAcoeff[3])*_indpndt_var+C_),2))),_ex1_2)
                };
    }
    else if(_odetype == odetype5)
    {
        const exset symClt = symbols(degAcoeff[1]+degAcoeff[3]+degAcoeff[5]+degAcoeff[7]);

        if(symClt.size()>1) // minimum two symbols are required for the following solutions with conditions, a0=8*a2^2/(27*a4) and a6=a4^2/(4*a2)
          {
               return {lst{pow((-8*degAcoeff[3]*pow(tanh(sqrt(-(_ex1/3)*degAcoeff[3])*_indpndt_var+C_),2)/(3*degAcoeff[5]*(3+pow(tanh(sqrt(-(_ex1/3)*degAcoeff[3])*_indpndt_var+C_),2)))),_ex1_2),
                       degAcoeff[1]==8*pow(degAcoeff[3],2)/(27*degAcoeff[5]),degAcoeff[7]==pow(degAcoeff[5],2)/(4*degAcoeff[3])},
                       lst{pow((-8*degAcoeff[3]*pow(coth(sqrt(-(_ex1/3)*degAcoeff[3])*_indpndt_var+C_),2)/(3*degAcoeff[5]*(3+pow(coth(sqrt(-(_ex1/3)*degAcoeff[3])*_indpndt_var+C_),2)))),_ex1_2),
                       degAcoeff[1]==8*pow(degAcoeff[3],2)/(27*degAcoeff[5]),degAcoeff[7]==pow(degAcoeff[5],2)/(4*degAcoeff[3])}
                      };
          }

    }
    else if(_odetype == odetype6)
    {
        return {-degAcoeff[3]*degAcoeff[4]*pow(sech((_ex1_2)*sqrt(degAcoeff[3])*_indpndt_var),2)/(pow(degAcoeff[4],2)-degAcoeff[3]*degAcoeff[5]*pow((1-tanh((_ex1_2)*sqrt(degAcoeff[3])*_indpndt_var)),2)),
               2*degAcoeff[3]*sech(sqrt(degAcoeff[3])*_indpndt_var)/(sqrt(-4*degAcoeff[3]*degAcoeff[5]+pow(degAcoeff[4],2))-degAcoeff[4]*sech(sqrt(degAcoeff[3])*_indpndt_var))};
    }
    else if(_odetype == odetype7)
    {
        const ex A0A1 = simplify(degAcoeff[1]/degAcoeff[2]);
        return {pow((C_ + _indpndt_var),2)*degAcoeff[2]/(_ex2*_ex2) - A0A1};
    }
    else if (_odetype == jacobiElip024)
    {
        const ex discrimi = simplify(-4*degAcoeff[1]*degAcoeff[5] + pow(degAcoeff[3],2));

        return {sqrt(-degAcoeff[3]*degAcoeff[5]*(pow(degAcoeff[3],2) - sqrt(pow(degAcoeff[3],2)*(discrimi))))*JacobiSN(sqrt(degAcoeff[3]*(-pow(degAcoeff[3],2) - sqrt(pow(degAcoeff[3],2)*(discrimi))))*_indpndt_var/(sqrt(_ex2)*degAcoeff[3]), sqrt(degAcoeff[1]*degAcoeff[5]*(-2*degAcoeff[1]*degAcoeff[5] + pow(degAcoeff[3],2) - sqrt(pow(degAcoeff[3],2)*(discrimi))))/(sqrt(_ex2)*degAcoeff[1]*degAcoeff[5]))/(sqrt(_ex2)*degAcoeff[3]*degAcoeff[5]),
                sqrt(degAcoeff[3]*degAcoeff[5]*(-pow(degAcoeff[3],2) + sqrt(pow(degAcoeff[3],2)*(discrimi))))*JacobiCN(sqrt(-degAcoeff[3]*sqrt(pow(degAcoeff[3],2)*(discrimi)))*_indpndt_var/degAcoeff[3], sqrt(-(discrimi)*(sqrt(pow(degAcoeff[3],2)*(discrimi)) - (discrimi)))/(sqrt(_ex2)*(discrimi)))/(sqrt(_ex2)*degAcoeff[3]*degAcoeff[5]),
                sqrt(-degAcoeff[3]*degAcoeff[5]*(pow(degAcoeff[3],2) + sqrt(pow(degAcoeff[3],2)*(discrimi))))*JacobiDN(sqrt(degAcoeff[3]*(pow(degAcoeff[3],2) + sqrt(pow(degAcoeff[3],2)*(discrimi))))*_indpndt_var/(sqrt(_ex2)*degAcoeff[3]), sqrt(degAcoeff[1]*degAcoeff[5]*(sqrt(pow(degAcoeff[3],2)*(discrimi)) - (discrimi)))/(sqrt(_ex2)*degAcoeff[1]*degAcoeff[5]))/(sqrt(_ex2)*degAcoeff[3]*degAcoeff[5]),
                sqrt(degAcoeff[5]*(-degAcoeff[3] + sqrt(discrimi)))*JacobiNS(_ex1_2*sqrt(-2*degAcoeff[3] + 2*sqrt(discrimi))*_indpndt_var, sqrt(degAcoeff[1]*degAcoeff[5]*(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[1]*degAcoeff[5] + pow(degAcoeff[3],2)))/(sqrt(_ex2)*degAcoeff[1]*degAcoeff[5]))/(sqrt(_ex2)*degAcoeff[5]),
                sqrt(-degAcoeff[3] - sqrt(discrimi))*JacobiNC(sqrt(-sqrt(discrimi))*_indpndt_var, sqrt(-2*degAcoeff[1]*degAcoeff[5])/sqrt(degAcoeff[3]*sqrt(discrimi) + discrimi))/sqrt(2*degAcoeff[5]),
                sqrt(-degAcoeff[3] + sqrt(discrimi))*JacobiND(sqrt(-2*degAcoeff[1]*degAcoeff[5])*_indpndt_var/sqrt(-degAcoeff[3] + sqrt(discrimi)), sqrt(degAcoeff[3]*sqrt(discrimi) - discrimi)/sqrt(2*degAcoeff[1]*degAcoeff[5]))/sqrt(2*degAcoeff[5]),
                sqrt(2*degAcoeff[1])*JacobiSC(-sqrt(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[1]*degAcoeff[5] + pow(degAcoeff[3],2))*_indpndt_var/sqrt(degAcoeff[3] + sqrt(discrimi)), sqrt(degAcoeff[3]*sqrt(discrimi) + discrimi)/sqrt(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[1]*degAcoeff[5] + pow(degAcoeff[3],2)))/sqrt(degAcoeff[3] + sqrt(discrimi)),
                sqrt(2*degAcoeff[1])*sqrt(2*degAcoeff[1]*degAcoeff[5] - pow(degAcoeff[3],2) - degAcoeff[3]*sqrt(discrimi))*JacobiSD(sqrt(degAcoeff[3]*sqrt(discrimi) + discrimi)*_indpndt_var/sqrt(degAcoeff[3] + sqrt(discrimi)), sqrt(2*degAcoeff[1]*degAcoeff[5] - pow(degAcoeff[3],2) - degAcoeff[3]*sqrt(discrimi))/sqrt(-discrimi - degAcoeff[3]*sqrt(discrimi)))/(sqrt(degAcoeff[3] + sqrt(discrimi))*sqrt(-discrimi - degAcoeff[3]*sqrt(discrimi))),
                sqrt(2*degAcoeff[1])*JacobiCS(sqrt(2*degAcoeff[1]*degAcoeff[5])*_indpndt_var/sqrt(degAcoeff[3] + sqrt(discrimi)), sqrt(-degAcoeff[3]*sqrt(discrimi) - discrimi)/sqrt(2*degAcoeff[1]*degAcoeff[5]))/sqrt(degAcoeff[3] + sqrt(discrimi)),
                sqrt(degAcoeff[3]*sqrt(discrimi) - discrimi)*JacobiDS(sqrt(degAcoeff[3]*sqrt(discrimi) - discrimi)*_indpndt_var/sqrt(degAcoeff[3] - sqrt(discrimi)), sqrt(2*degAcoeff[1]*degAcoeff[5])/sqrt(degAcoeff[3]*sqrt(discrimi) - discrimi))/sqrt(degAcoeff[3]*degAcoeff[5] - degAcoeff[5]*sqrt(discrimi))
                };

    }
    else if (_odetype==jacobiElip123)
    {
        const ex discrimi = simplify(-4*degAcoeff[2]*degAcoeff[4]+pow(degAcoeff[3],2));

        return {_ex1_2*(-degAcoeff[3] + sqrt(discrimi))*pow(JacobiSN(sqrt(degAcoeff[2]*degAcoeff[4])*_indpndt_var/(sqrt(_ex2)*sqrt(-degAcoeff[3] + sqrt(discrimi))), sqrt(pow(degAcoeff[3],2) - degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[2]*degAcoeff[4])/(sqrt(_ex2)*sqrt(degAcoeff[2]*degAcoeff[4]))),2)/degAcoeff[4],
                -2*degAcoeff[2]*pow(JacobiCN(-pow((discrimi),(_ex1_2*_ex1_2))/2*_indpndt_var, sqrt(degAcoeff[3] + sqrt(discrimi))/(sqrt(_ex2)*pow((discrimi),(_ex1_2*_ex1_2)))),2)/(degAcoeff[3] - sqrt(discrimi)),
                -2*degAcoeff[2]*pow(JacobiDN(-sqrt(degAcoeff[2]*degAcoeff[4])*_indpndt_var/(sqrt(_ex2)*sqrt(degAcoeff[3] - sqrt(discrimi))), sqrt(4*degAcoeff[2]*degAcoeff[4] + degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))/sqrt(2*degAcoeff[2]*degAcoeff[4])),2)/(degAcoeff[3] - sqrt(discrimi)),
                -2*degAcoeff[2]*pow(JacobiNS(sqrt(-degAcoeff[3] + sqrt(discrimi))*_indpndt_var/(2*sqrt(_ex2)), sqrt(degAcoeff[3] + sqrt(discrimi))/sqrt(degAcoeff[3] - sqrt(discrimi))),2)/(degAcoeff[3] + sqrt(discrimi)),
                -2*degAcoeff[2]*pow(JacobiNC(-pow((discrimi),(_ex1_2*_ex1_2))/2*_indpndt_var, sqrt(degAcoeff[3] + sqrt(discrimi))/(sqrt(_ex2)*pow((discrimi),(_ex1_2*_ex1_2)))),2)/(degAcoeff[3] + sqrt(discrimi)),
                -2*degAcoeff[2]*pow(JacobiND(-sqrt(degAcoeff[3] + sqrt(discrimi))*_indpndt_var/(2*sqrt(_ex2)), -sqrt(_ex2)*pow((discrimi),(_ex1_2*_ex1_2))/sqrt(degAcoeff[3] + sqrt(discrimi))),2)/(degAcoeff[3] + sqrt(discrimi)),
                2*degAcoeff[2]*pow(JacobiSC(sqrt(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))*_indpndt_var/(2*sqrt(degAcoeff[3] + sqrt(discrimi))), sqrt(degAcoeff[3]*sqrt(discrimi) - 4*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))/sqrt(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))),2)/(degAcoeff[3] + sqrt(discrimi)),
                2*degAcoeff[2]*(degAcoeff[3]*sqrt(discrimi) - 2*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))*pow(JacobiSD(sqrt(degAcoeff[3]*sqrt(discrimi) - 4*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))*_indpndt_var/(2*sqrt(degAcoeff[3] + sqrt(discrimi))), sqrt(2*degAcoeff[2]*degAcoeff[4] - degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))/sqrt(4*degAcoeff[2]*degAcoeff[4] - degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))),2)/((degAcoeff[3] + sqrt(discrimi))*(degAcoeff[3]*sqrt(discrimi) - 4*degAcoeff[2]*degAcoeff[4] + pow(degAcoeff[3],2))),
                2*degAcoeff[2]*pow(JacobiCS(sqrt(2*degAcoeff[2]*degAcoeff[4])*_indpndt_var/(2*sqrt(degAcoeff[3] - sqrt(discrimi))), sqrt(4*degAcoeff[2]*degAcoeff[4] + degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))/sqrt(2*degAcoeff[2]*degAcoeff[4])),2)/(degAcoeff[3] - sqrt(discrimi)),
                (4*degAcoeff[2]*degAcoeff[4] - degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))*pow(JacobiDS(sqrt(4*degAcoeff[2]*degAcoeff[4] - degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))*_indpndt_var/(2*sqrt(degAcoeff[3] + sqrt(discrimi))), sqrt(2*degAcoeff[2]*degAcoeff[4])/sqrt(4*degAcoeff[2]*degAcoeff[4] - degAcoeff[3]*sqrt(discrimi) - pow(degAcoeff[3],2))),2)/((degAcoeff[3] + sqrt(discrimi))*degAcoeff[4])

                };

    }

    return {};
}

/**checking independent vars in diff. equ. **/
ex diffIndpndtPrsntTest::operator()(const ex& _e)
{
    if(is_ex_the_function(_e, Diff))
    {
        return _e;
    }
    else if(is_a<symbol>(_e))
    {
        for(auto it = indpndtVarLst.begin(); it != indpndtVarLst.end(); it++)
        {
            if(_e==*it)
                isIndpndtPrsnt = true;
        }
        return _e;
    }

    return (_e.map(*this));
}

/// depended variable (_y) is replaced by _y^_n.
/// negative signs in power are turned into positive signs.
/// summing powers of a term and the summed power are collected in n_pow_clt.
/// this is repeated for other terms.
int find_n_power::n_pow(const ex& expr_, const ex& dpndt_var)
{
    ex _pairpop;
    exmap repls;
    int ret;

    // y and $*y
    if((expr_) == dpndt_var)
    {
        if(multindic == 0)
            n_pow_clt.push_back(_n);
        else
        {
            if(!n_pow_clt.empty())
            {
                _pairpop = n_pow_clt.back();
                n_pow_clt.pop_back();
                n_pow_clt.push_back(_pairpop + _n);
            }
        }
        return 0;
    }
    repls.clear();

    // y^2 and $*y^2
    (expr_).match(pow(dpndt_var, wild(0)), repls);
    if(!repls.empty())
    {
        if(multindic == 0)
            n_pow_clt.push_back(_n*repls[wild(0)]);
        else
        {
            if(!n_pow_clt.empty())
            {
                _pairpop = n_pow_clt.back();
                n_pow_clt.pop_back();
                n_pow_clt.push_back( _pairpop + _n*repls[wild(0)]);
            }
        }

        return 0;
    }
    repls.clear();

    // diff(y, x, 1) and $*diff(y, x, 1)
    (expr_).match(Diff(dpndt_var, wild(0), wild(1)), repls);

    if(!repls.empty())
    {
        if(multindic == 0)
            n_pow_clt.push_back(_n+repls[wild(1)]);
        else
        {
            if(!n_pow_clt.empty())
            {
                _pairpop = n_pow_clt.back();
                n_pow_clt.pop_back();
                n_pow_clt.push_back( _pairpop + _n+repls[wild(1)]);
            }
        }

        return 0;
    }
    repls.clear();

    // diff(y, x, 1)^2 and $*diff(y, x, 1)^2
    (expr_).match(pow(Diff(dpndt_var, wild(0), wild(1)), wild(2)), repls);
    if(!repls.empty())
    {
        if(multindic == 0)
            n_pow_clt.push_back(repls[wild(2)]*(_n+repls[wild(1)]));
        else
        {
            if(!n_pow_clt.empty())
            {
                _pairpop = n_pow_clt.back();
                n_pow_clt.pop_back();
                n_pow_clt.push_back( _pairpop + repls[wild(2)]*(_n+repls[wild(1)]));
            }
        }

        return 0;
    }
    repls.clear();

    if(is_a<add>(expr_))
    {
        for(size_t i=0; i<expr_.nops(); i++)
        {
            ret = this->n_pow(expr_.op(i), dpndt_var);

            if(ret == -1 && is_a<mul>(expr_.op(i)))
            {
               for(size_t j=0; j<(expr_.op(i)).nops(); j++)
               {
                  ret = this->n_pow(expr_.op(i).op(j), dpndt_var);
                  if(ret == 0)
                    multindic = 1;
               }
               multindic = 0;
            }
        }
    }

    return -1;
}

/// finding order of ode
/////////////////////////////////////////////////////////////////////
ex find_ode_order::operator()(const ex& _e)
{
   if(is_ex_the_function(_e, Diff))
   {
        if(!(_e.op(2)).info(info_flags::integer))
        {
            dorat dorat; // avoiding double form
            dorat.set();
            if(!dorat(_e.op(2)).info(info_flags::integer))
                return _FAIL;
            order_clt.push_back(ex_to<numeric>(dorat(_e.op(2))).to_int());
        }
        else
            order_clt.push_back(ex_to<numeric>(_e.op(2)).to_int());

        return _e.map(*this);
   }
   return _e.map(*this);
}
///////////////////////////////////////////////////////////////////////

int order(const ex& _expr)
{
    find_ode_order find_ode_order;
    find_ode_order(_expr);
    if(!find_ode_order.order_clt.empty())
        return *max_element(find_ode_order.order_clt.begin(), find_ode_order.order_clt.end());
    else
        return 0;
}
///////////////////////////////////////////////////////////////////////

ex find_indpndt_var::operator()(const ex& _e)
{
    if(is_ex_the_function(_e, Diff))
        indpndt_var.insert(_e.op(1));
    return _e.map(*this);
}

//////////////////////////////////////////////////////////////////////////

ex twf::operator()(const ex& _e)
{
    if(is_ex_the_function(_e, Diff) && _e.op(1) != xi)
    {
        return Diff(_e.op(0), xi, _e.op(2))*pow(indpndt_vars_wt_indx[_e.op(1)], _e.op(2));
    }
    else
    {
        return _e.map(*this);
    }
}

//////////////////////////////////////////////////////////////////////////////

int desolve(const ex& diffeq, const lst& dpndt_vars, const int& method, bool test)
{

    constraints.remove_all();
    solutionClt.clear();

    stringstream solutions, strStrm;
    ex temdiffeq, remainingDiffpart = _ex0;
    // collecting variables from twc, phase, paraInDiffSolve for which system of equations to be solved.
    lst variables;


    vector<string> symbolClt, protectedSymbols;
    stringstream symbolStr;

    /*list of protected symbols*/
    protectedSymbols.push_back("N");
    protectedSymbols.push_back("F");
    protectedSymbols.push_back("Fd");
    protectedSymbols.push_back("Xun");
    protectedSymbols.push_back("Yun");
    protectedSymbols.push_back("hun");
    protectedSymbols.push_back("gun");
    protectedSymbols.push_back("a0Deg");
    protectedSymbols.push_back("gDeg");
    protectedSymbols.push_back("a1Deg");
    protectedSymbols.push_back("Const");
    protectedSymbols.push_back("F_");
    protectedSymbols.push_back("Fd_");
    protectedSymbols.push_back("X_");
    protectedSymbols.push_back("Y_");
    protectedSymbols.push_back("h_");
    protectedSymbols.push_back("g_");
    protectedSymbols.push_back("C_");
    for(unsigned i=0;i<10;i++)
    {

        symbolStr<<"a_"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");
        symbolStr<<"b_"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");
        symbolStr<<"g_"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");

        symbolStr<<"a"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");
        symbolStr<<"b"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");
        symbolStr<<"g"<<i;
        protectedSymbols.push_back(symbolStr.str());
        symbolStr.str("");
    }
    for(unsigned i=0;i<10;i++)
    {
        for(unsigned j=0;j<10;j++)
        {

            symbolStr<<"a_"<<i<<j;
            protectedSymbols.push_back(symbolStr.str());
            symbolStr.str("");
            symbolStr<<"b_"<<i<<j;
            protectedSymbols.push_back(symbolStr.str());
            symbolStr.str("");

            symbolStr<<"a"<<i<<j;
            protectedSymbols.push_back(symbolStr.str());
            symbolStr.str("");
            symbolStr<<"b"<<i<<j;
            protectedSymbols.push_back(symbolStr.str());
            symbolStr.str("");
        }
    }


   symbol_finder.clear();
   symbol_finder(diffeq);
   for(auto itr=symbol_finder.symbols.begin();itr!=symbol_finder.symbols.end();itr++)
   {
       symbolStr<<*itr;
       symbolClt.push_back(symbolStr.str());
       symbolStr.str("");
   }
   symbol_finder.clear();

   for(auto itr=protectedSymbols.begin();itr!=protectedSymbols.end();itr++)
   {
       if(std::find(symbolClt.begin(),symbolClt.end(),*itr)!=symbolClt.end())
       {
           cout<<"Evaluation stop: "<<*std::find(symbolClt.begin(),symbolClt.end(),*itr)<<" symbol is protected"<<endl;

           solutions <<"Evaluation stop: "<<*std::find(symbolClt.begin(),symbolClt.end(),*itr)<<" symbol is protected"<<endl;
           writetofile(solutions, dpndt_vars.op(0));

           #ifdef GiNaCDE_gui
           stringstream temstr;
           temstr <<"Evaluation stop: "<<*std::find(symbolClt.begin(),symbolClt.end(),*itr)<<" symbol is protected.";
           gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, &temstr.str()[0]);
           //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
           #endif // GiNaCDE_gui
           return -1;
       }
   }


    if(method != FIM && nops(degAcoeff) == 0)
    {
        cout <<"Evaluation stop: Please initialize 'degAcoeff';"<<endl;
        solutions << "Evaluation stop: Please initialize 'degAcoeff';"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Please initialize 'degAcoeff';");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }

    if(method != FIM && nops(degAcoeff) != 0 && !is_a<numeric>(degAcoeff[0]))
    {
        cout <<"Evaluation stop: Please provide highest positive integer delta of 1st order NLODE in 'degAcoeff';"<<endl;
        solutions << "Evaluation stop: Please provide highest positive integer delta of 1st order NLODE in 'degAcoeff';"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Please provide highest positive integer delta of 1st order NLODE in 'degAcoeff';");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }


    if(method != FIM && !degAcoeff[0].info(info_flags::positive))
    {
        cout <<"Evaluation stop: delta of 1st order NLODE in 'degAcoeff' should be positive integer;"<<endl;
        solutions << "Evaluation stop: delta of 1st order NLODE in 'degAcoeff' should be positive integer;"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: delta of 1st order NLODE in 'degAcoeff' should be positive integer;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }

    if(method != FIM && ex_to<numeric>(degAcoeff[0]).to_int()!=(int)(nops(degAcoeff)-2))
    {
        cout <<"Evaluation stop: Please provide all coefficients of 1st order NLODE in 'degAcoeff';"<<endl;
        solutions << "Evaluation stop: Please provide all coefficients of 1st order NLODE in 'degAcoeff';"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Please provide all coefficients of 1st order NLODE in 'degAcoeff';");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }

    if(!positivePart && !negativePart)
    {
        cout <<"Evaluation stop: Minimum one part (positivePart or negativePart) should be \"yes\";"<<endl;
        solutions << "Evaluation stop: Minimum one part (positivePart or negativePart) should be \"yes\";"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        return -1;
    }

    if(method == FIM && NValue!=1 && NValue!=2 && NValue!=0)
    {
        cout <<"Evaluation stop: For FIM method, value of N is either 1 or 2';"<<endl;
        solutions << "Evaluation stop: For FIM method, value of N is either 1 or 2';"<<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: For FIM method, value of N is eighter 1 or 2';");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }

    //temdiffeq = Simplify(evaluate(Simplify(diffeq)));

    temdiffeq = Simplify(expand(Simplify(diffeq)));

    if(temdiffeq == _FAIL)
    {
        cout << "Evaluation stop: unsupported diff. Equ.;" << endl;
        solutions << "Evaluation stop: unsupported diff. Equ.;" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported diff. Equ.;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui

        return -1;
    }


    if(!temdiffeq.has(Diff(wild(0),wild(1),wild(2))))
    {
        cout << "Evaluation stop: it is not a diff. Equ.;" << endl;
        solutions << "Evaluation stop: it is not a diff. Equ.;" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: It is not a diff. Equ.;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui

        return -1;
    }

    // handling denominator of ode //
    lst diffDenomSolu = {};
    if(denom(temdiffeq) != _ex1 && Simplify(expand(denom(temdiffeq))).has(dpndt_vars[0]))
    {
        const exsetlst diffDenomSoluC = solve({Simplify(expand(denom(temdiffeq)))}, {dpndt_vars[0]});

        if(!diffDenomSoluC.empty())
        {
            for(exsetlst::const_iterator it1 = diffDenomSoluC.begin(); it1 != diffDenomSoluC.end(); it1++)
            {
                diffDenomSolu.append((*it1)[0].rhs());

            }
        }

    }

    temdiffeq = numer(temdiffeq);






    if( method == FIM )
        solutions << "                              ==========First Integral Method===========                         " << endl;
    else if( method == F_expansion )
        solutions << "                              ==========F-Expansion Method===========                         " << endl;
    else if( method == mF_expansion )
        solutions << "                              ==========Modified F-Expansion Method===========                         " << endl;

    if(output==maple )
        solutions <<  "Equations are written in MAPLE language." << endl;
    else if(output==mathematica)
        solutions <<  "Equations are written in MATHEMATICA language." << endl;
    else
        solutions <<  "Equations are written in GiNaC language." << endl;

    const string out = outstr("-", 100);
    solutions << out << endl << endl;

    xi=reader("xi"); // traveling wave coordinate
    bool phasepart = false;

    ex exTem;

    find_indpndt_var find_indpndt_var;
    find_indpndt_var(temdiffeq);
    exset indpndt_vars = find_indpndt_var.indpndt_var;

    /*checking independent vars in diff. equ.*/
    diffIndpndtPrsntTest diffIndpndtPrsntTest(false,indpndt_vars);
    diffIndpndtPrsntTest(temdiffeq);
    if(diffIndpndtPrsntTest.isIndpndtPrsnt)
    {
        cout << "Evaluation stop: independent vars present in diff. Equ.;" << endl;
        solutions << "Evaluation stop: independent vars present in diff. Equ.;" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: independent vars present in diff. Equ.;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }

    // arranging independent variables //
    exvector indpndt_varsVec;
    for( auto it = var_depend[dpndt_vars.op(0)].begin(); it != var_depend[dpndt_vars.op(0)].end(); it++ )
    {
        if( find(indpndt_vars.begin(), indpndt_vars.end(), *it) != indpndt_vars.end() )
           {if( find(indpndt_varsVec.begin(), indpndt_varsVec.end(), *it) == indpndt_varsVec.end() ) indpndt_varsVec.push_back( *it );}
    }

    if(indpndt_vars.size()!=indpndt_varsVec.size())
    {
        cout << "Evaluation stop: Please provide dependency of all independent variables" <<endl;
        solutions << "Evaluation stop: Please provide dependency of all independent variables" <<endl;
        writetofile(solutions, dpndt_vars.op(0));
        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Please provide dependency of all independent variables");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }



    ex U=reader("U");

    for(auto it = indpndt_varsVec.begin(); it != indpndt_varsVec.end(); it++ )
    {
        depend(U, {*it});
    }


    cout<<endl<<endl;
    solutions << "Input equation is: ";
    cout << "Input equation is: ";

    if(output == maple)
    {
        solutions << diffformchange(diffeq, dpndt_vars, indpndt_vars) << " = 0;" << endl;
        cout << diffformchange(diffeq, dpndt_vars, indpndt_vars) << " = 0;" << endl;
    }
    else if(output == mathematica)
    {
        solutions <<  gmathematica(diffformchange(diffeq, dpndt_vars, indpndt_vars)) << " = 0;" << endl;
        cout <<  gmathematica(diffformchange(diffeq, dpndt_vars, indpndt_vars)) << " = 0;" << endl;
    }
    else
    {
        solutions << diffeq << " = 0;" << endl;
        cout << diffeq << " = 0;" << endl;
    }



    // converting pde to twf(travelling wave form)
    if( indpndt_vars.size() > 1 && nops(twcPhase)!=2)
    {
        cout << "Evaluation stop: please use 'twcPhase' properly;" << endl;
        solutions << "Evaluation stop: please use 'twcPhase' properly;" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: please use 'twcPhase' properly;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }
    lst twc, phase;

    ex tw_coordi = _ex0, tw_coordiPhase = _ex0;

    string symname;
    //conjuFree conjuFree;
    replaceI replaceI;
    ex replaceIex;

    if(indpndt_vars.size() > 1)
    {
#ifndef GiNaCDE_gui
        twc=ex_to<lst>(twcPhase[0]), phase=ex_to<lst>(twcPhase[1]);
#endif //GiNaCDE_gui

        for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
        {
            if((temdiffeq.subs(Diff(wild(0), wild(1), wild(2)) == xi)).has(*it))
            {
                cout << "Evaluation stop: unsupported diff. Equ.;" << endl;
                solutions << "Evaluation stop: unsupported diff. Equ.;" << endl;
                writetofile(solutions, dpndt_vars.op(0));

                #ifdef GiNaCDE_gui
                gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported diff. Equ.;");
                //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                #endif // GiNaCDE_gui
                return -1;
            }
        }

        //Creating wave vector(k) & velocity(v)

        string tem1;

        #ifdef GiNaCDE_gui
        string tem;

        GtkWidget *table = gtk_grid_new();

        entry.clear();
        entrylbl.clear();
        entryText.clear();
        stringstream exstring;

        GtkWidget *twPartlbl = gtk_label_new("Constants in traveling wave part:");
        gtk_grid_attach(GTK_GRID(table),twPartlbl,0,0,2,1);
        for(unsigned i=0; i<indpndt_varsVec.size(); i++)
        {
            exstring.str("");
            exstring << "Constant with " << indpndt_varsVec[i]  << " coordinate: ";
            entrylbl.push_back(gtk_label_new(&(exstring.str()[0])));
            entry.push_back(gtk_entry_new());
            symname = "k_" + to_string(i);
            gtk_entry_set_text(GTK_ENTRY(entry[i]),&symname[0]);
            gtk_grid_attach(GTK_GRID(table),entrylbl[i],0,1+i,1,1);
            gtk_grid_attach(GTK_GRID(table),entry[i],1,1+i,1,1);
        }


        GtkWidget *dialog=gtk_dialog_new();
        gtk_window_set_decorated(GTK_WINDOW(dialog),FALSE);
        gtk_window_set_title(GTK_WINDOW(dialog),"Constants in traveling wave part:");
        gtk_window_set_modal( GTK_WINDOW(dialog), TRUE );
        gtk_window_set_transient_for( GTK_WINDOW(dialog), GTK_WINDOW(window));

        GtkAccelGroup *accel_group=gtk_accel_group_new();
        GtkWidget *okbutton=gtk_button_new_with_label("Ok");
        GtkWidget *cancelbutton=gtk_button_new_with_label("Cancel");
        gtk_window_add_accel_group(GTK_WINDOW(dialog),accel_group);
        gtk_widget_add_accelerator(okbutton,"clicked",accel_group,GDK_KEY_Return,(GdkModifierType)0, (GtkAccelFlags)0);
        GtkWidget *hbox=gtk_box_new(GTK_ORIENTATION_HORIZONTAL,10);
        gtk_box_pack_start(GTK_BOX(hbox),okbutton,FALSE,FALSE,0);
        gtk_box_pack_start(GTK_BOX(hbox),cancelbutton,FALSE,FALSE,0);
        gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog)))),table);
        //gtk_dialog_add_action_widget(GTK_DIALOG(dialog), hbox, GTK_RESPONSE_NONE);
        gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog)))),hbox);


        g_signal_connect(okbutton,"clicked",G_CALLBACK(twc_okbutton_clicked),(gpointer)dialog);
        g_signal_connect(cancelbutton,"clicked",G_CALLBACK(cancelbutton_clicked),(gpointer)dialog);
        for(unsigned i=0; i<entry.size(); i++)
            g_signal_connect(entry[i], "key-press-event", G_CALLBACK(key_press_cb),(gpointer)entrylbl[i]);

        gtk_widget_show_all(dialog);
        gtk_dialog_run(GTK_DIALOG(dialog));

        if(cancel)
        {
            cancel = false;
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
            return 0;
        }

        if(!entryText.empty())
        {
            for(unsigned i = 0; i < entryText.size(); i++)
            {
                twc.append(entryText[i]);
            }
        }
        #endif // GiNaCDE_gui

        if(nops(twc)==0)
        {
            cout << "Evaluation stop: please provide constants in traveling wave coordinate(twc);" << endl;
            solutions << "Evaluation stop: please provide constants in traveling wave coordinate(twc);" << endl;
            writetofile(solutions, dpndt_vars.op(0));

            #ifdef GiNaCDE_gui
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: please provide constants in traveling wave coordinates(twc);");
            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
            #endif // GiNaCDE_gui
            return -1;
        }

        if(nops(twc)!=indpndt_varsVec.size())
        {
            cout << "Evaluation stop: please provide constants for each independent variables in twc part of 'twcPhase';" << endl;
            solutions << "Evaluation stop: please provide constants for each independent variables in twc part of 'twcPhase';" << endl;
            writetofile(solutions, dpndt_vars.op(0));

            #ifdef GiNaCDE_gui
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: please provide constants for each independent variables in twc part of 'twcPhase';");
            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
            #endif // GiNaCDE_gui
            return -1;
        }

        std::map<ex, ex, ex_is_less> indpndt_vars_wt_indx;
        size_t i = 0;

        for(auto it = indpndt_varsVec.begin(); it != indpndt_varsVec.end(); it++ )
        {
            indpndt_vars_wt_indx[*it] = twc[i];
            i = i + 1;
        }

        replaceIex = replaceI(temdiffeq);
        const ex conjuFreeex = conjuFreee(temdiffeq);
        if(replaceIex != temdiffeq || conjuFreeex != temdiffeq) //complex NLPDE
        {
            #ifndef GiNaCDE_gui
            if(nops(phase)==0)
            {
                cout <<  "Evaluation stop: please provide phase part for complex NLPDE;"  << endl;
                solutions << "Evaluation stop: please provide phase part for complex NLPDE;" << endl;
                writetofile(solutions, dpndt_vars.op(0));
                phasepart=false;
                return -1;
            }
            else
            {
                if(nops(phase)!=nops(twc))
                {
                    cout <<  "Evaluation stop: please provide constants for each independent variables in phase part of 'twcPhase';"  << endl;
                    solutions << "Evaluation stop: please provide constants for each independent variables in phase part of 'twcPhase';" << endl;
                    writetofile(solutions, dpndt_vars.op(0));
                    return -1;
                }
                phasepart=true;
            }
            #endif // GiNaCDE_gui

            #ifdef GiNaCDE_gui
            string tem;

            table = gtk_grid_new();

            entry.clear();
            entrylbl.clear();
            entryText.clear();
            stringstream exstring;
            for(unsigned i=0; i<indpndt_varsVec.size(); i++)
            {
                exstring.str("");
                exstring << "Constant with " << indpndt_varsVec[i]  << " coordinate in phase part: ";
                entrylbl.push_back(gtk_label_new(&exstring.str()[0]));
                entry.push_back(gtk_entry_new());
                symname = "p_" + to_string(i);
                gtk_entry_set_text(GTK_ENTRY(entry[i]),&symname[0]);
                gtk_grid_attach(GTK_GRID(table),entrylbl[i],0,i,1,1);
                gtk_grid_attach(GTK_GRID(table),entry[i],1,i,1,1);
            }

            dialog=gtk_dialog_new();
            gtk_window_set_decorated(GTK_WINDOW(dialog),FALSE);
            gtk_window_set_title(GTK_WINDOW(dialog),"Phase part in diff. equ.:");
            gtk_window_set_modal( GTK_WINDOW(dialog), TRUE );
            gtk_window_set_transient_for( GTK_WINDOW(dialog), GTK_WINDOW(window));


            accel_group=gtk_accel_group_new();
            okbutton=gtk_button_new_with_label("Ok");
            cancelbutton=gtk_button_new_with_label("Cancel");
            gtk_window_add_accel_group(GTK_WINDOW(dialog),accel_group);
            gtk_widget_add_accelerator(okbutton,"clicked",accel_group,GDK_KEY_Return,(GdkModifierType)0, (GtkAccelFlags)0);
            hbox=gtk_box_new(GTK_ORIENTATION_HORIZONTAL,10);
            gtk_box_pack_start(GTK_BOX(hbox),okbutton,FALSE,FALSE,0);
            gtk_box_pack_start(GTK_BOX(hbox),cancelbutton,FALSE,FALSE,0);
            gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog)))),table);
            //gtk_dialog_add_action_widget(GTK_DIALOG(dialog), hbox, GTK_RESPONSE_NONE);
            gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog)))),hbox);


          g_signal_connect(okbutton,"clicked",G_CALLBACK(phasepart_okbutton_clicked),(gpointer)dialog);
          //g_signal_connect(phaseCombo,"changed",G_CALLBACK(on_changed__phaseCombo),(gpointer)phaseCombolbl);
          g_signal_connect(cancelbutton,"clicked",G_CALLBACK(cancelbutton_clicked),(gpointer)dialog);

          for(unsigned i=0; i<entry.size(); i++)
            g_signal_connect(entry[i], "key-press-event", G_CALLBACK(key_press_cb),(gpointer)entrylbl[i]);

          gtk_widget_show_all(dialog);
          int result = gtk_dialog_run(GTK_DIALOG(dialog));

          if(cancel)
          {
              cancel = false;
              gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
              return 0;
          }

          switch(result)
          {
             case GTK_RESPONSE_DELETE_EVENT:
                  gtk_widget_destroy(dialog);
                  phasepart = false;
                  break;
          }

          if(!entryText.empty())
          {
              phasepart = true;
              for(unsigned i = 0; i < entryText.size(); i++)
              {
                    phase.append(entryText[i]);
              }
          }

            #endif // GiNaCDE_gui
        }

        if(nops(twc)!=0) // collecting variables from twc and phase for which system of equations to be solved.
        {
            for(auto it = twc.begin(); it != twc.end(); it++)
            {
                exset tem = symbols( *it );
                for( auto it1 = tem.begin(); it1 != tem.end(); it1++ )
                {
                    variables.append(*it1);
                }
            }
        }

        if( nops(phase)!=0 )
        {
            for(auto it = phase.begin(); it != phase.end(); it++)
            {
                exset tem = symbols( *it );
                for( auto it1 = tem.begin(); it1 != tem.end(); it1++ )
                {
                    variables.append(*it1);
                }
            }
        }

        unsigned indx = 0;
        if( phasepart )
        {
            // Constructing traveling wave coordinates for phase part //
            for(exvector::const_iterator it = indpndt_varsVec.begin(); it != indpndt_varsVec.end(); it++)
            {
                tw_coordiPhase = tw_coordiPhase + (*it)*phase[indx];
                indx = indx + 1;

            }

            temdiffeq = Simplify(evaluate(Simplify(expand(conjuFreee(evaluate(temdiffeq.subs(dpndt_vars.op(0) == U*exp((I)*(tw_coordiPhase)))))/(exp(I*(tw_coordiPhase)))))));

            if( (temdiffeq).has(exp(wild(0))))
            {
                cout << "Evaluation stop: unsupported diff. Equ.; " <<temdiffeq<< endl;
                solutions << "Evaluation stop: unsupported diff. Equ.;" << endl;
                writetofile(solutions, dpndt_vars.op(0));

                #ifdef GiNaCDE_gui
                gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported diff. Equ.;");
                //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                #endif // GiNaCDE_gui
                return -1;
            }


        }
        else
        {
            temdiffeq = temdiffeq.subs(dpndt_vars.op(0) == U);
            //temdiffeq = Simplify(evaluate( temdiffeq ));
        }

        /// Constructing traveling wave coordinates ///
        indx = 0;
        for(exvector::const_iterator it = indpndt_varsVec.begin(); it != indpndt_varsVec.end(); it++)
        {
            tw_coordi = tw_coordi + (*it)*twc[indx];
            indx = indx + 1;

        }

        twf twf(indpndt_vars_wt_indx);
        ex xprev;
        do
        {
            xprev = temdiffeq;
            temdiffeq = twf(temdiffeq);
        } while(xprev != temdiffeq);

        indpndt_vars.clear();
        indpndt_vars.insert(xi);
        depend(U, {xi});

       // temdiffeq = evaluate(temdiffeq);
    }
    else // ODE
    {
        U = dpndt_vars.op(0);
        twcPhase = lst{lst{},lst{}};
        twc = {};
        phase = {};
    }


    polycheckFim odecheck(true,*indpndt_vars.begin());

    odecheck(conjuFreee(replaceI(temdiffeq)));
    if(!odecheck.polytype)
    {

        cout << "Evaluation stop: unsupported diff. Equ.;" << endl;
        solutions << "Evaluation stop: unsupported diff. Equ.;" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: unsupported diff. Equ.;");
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui
        return -1;
    }


//////////////////////////////////////////////////////////////////////
    bool integrateDone = false;
    if(phasepart)
    {
        replaceIex = replaceI(temdiffeq);

        if(replaceIex.is_polynomial(symb_) && degree(replaceIex,symb_)<=1)
        {
            ex realpart = coeff(replaceIex,symb_,0), // getting real part of diff. equa.
               imagepart = coeff(replaceIex,symb_,1); // getting imaginary part of diff. equ.

            solutions << "Real part of Diff. Equ.: " <<diffformchange(realpart, (lst){U}, indpndt_vars) <<" = 0;"<< endl;
            solutions << "Imaginary part of Diff. Equ.: " <<diffformchange(imagepart, (lst){U}, indpndt_vars) <<" = 0;"<< endl;

            ex odeprev,temex;
            temex=realpart;
            bool cond = true;
            if(realpart.has(U)) // doing integration of real part.
            {
                do
                {
                    odeprev = integrate(realpart, *indpndt_vars.begin());

                    if(!odeprev.has(Integrate(wild(0), wild(1), wild(2))))
                    {

                        realpart = odeprev;
                    }
                    else
                        cond = false;

                }while(cond);
                if(realpart!=temex)
                {
                    solutions << "Real part of Diff. Equ. is integrable.\n"\
                                "After integration, Real part: " <<diffformchange(realpart, (lst){U}, indpndt_vars) <<" = 0;"<<endl;
                }
            }

            cond = true;
            temex = imagepart;
            if(imagepart.has(U))  // doing integration of imaginary part.
            {
                do
                {
                    odeprev = integrate(imagepart, *indpndt_vars.begin());

                    if(!odeprev.has(Integrate(wild(0), wild(1), wild(2))))
                    {

                        imagepart = odeprev;
                    }
                    else
                        cond = false;

                }while(cond);
                if(imagepart!=temex)
                {
                    solutions << "Imaginary part of Diff. Equ. is integrable.\n"\
                    "After integration, Imaginary part: " <<diffformchange(imagepart, (lst){U}, indpndt_vars) <<" = 0;"<< endl;
                }
            }

            if(imagepart == _ex0 || realpart == _ex0)
            {

                temdiffeq=realpart+I*imagepart;
                solutions << "The Diff. Equ. becomes: " << diffformchange(temdiffeq, (lst){U}, indpndt_vars) <<" = 0;" << endl;
            }
            else
            {
                exset temSymb = symbols(imagepart); // checking present of symbols in imagepart
                if(temSymb.empty())
                {
                    cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                    solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                    writetofile(solutions, dpndt_vars.op(0));

#ifdef GiNaCDE_gui
                    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
                    //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
#endif // GiNaCDE_gui

                    return -1;
                }

                cond = false;
                temex = collect_common_factors(imagepart);
                ex temexClt=_ex1;

                if (is_a<mul>(temex))
                {
                    for (size_t i = 0; i < nops(temex);i++)
                    {
                        if (!temex.op(i).has(U) && is_a<add>(temex.op(i)))
                        {
                            cond = true;
                            temexClt = temexClt*temex.op(i);
                        }
                    }
                }
                else if(!temex.has(U))
                {
                    cond = true;
                }
                temex = temexClt;

                if(cond)
                {
                    temSymb = symbols(temex);
                    vector<ex> varSlct; // selecting parameter to get condition from imaginary part
                    if(temSymb.size()>1)
                    {
                        std::vector<std::pair<ex, int>> temSymbMapWtDeg;
                        for (auto it = temSymb.begin();it != temSymb.end();it++)
                        {
                            if(is_polynomial(temex,*it))
                                temSymbMapWtDeg.push_back(make_pair(*it,degree(temex,*it)));
                            else
                                varSlct.push_back(*it);
                        }
                        if(temSymbMapWtDeg.size()>1)
                        {
                            sort(temSymbMapWtDeg.begin(), temSymbMapWtDeg.end(), [=](std::pair<ex, int>& a, std::pair<ex, int>& b){
                                return a.second < b.second;});
                            varSlct.clear();
                            varSlct.push_back(temSymbMapWtDeg.begin()->first);
                        }
                    }
                    else if(!temSymb.empty())
                        varSlct.push_back(*temSymb.begin());

                    exsetlst imagepartSolu;

                    if(!varSlct.empty())
                        imagepartSolu = solve({temex},{*varSlct.begin()});
                    if(!imagepartSolu.empty())
                    {
                        temdiffeq = realpart;
                        remainingDiffpart = imagepart;
                        solutions << "We derive solutions of real part of Diff. Equ. with the condition:"
                                     "\n "<<(*imagepartSolu.begin()).op(0).lhs()<<" = "<<(*imagepartSolu.begin()).op(0).rhs()<<";"<<endl;

                        constraints.append((*imagepartSolu.begin()).op(0).lhs()==(*imagepartSolu.begin()).op(0).rhs());

                        ex twcex = twc,  phaseex = phase; // substitute (*imagepartSolu.begin()).op(0).lhs() = (*imagepartSolu.begin()).op(0).rhs() in twc, phase, temdiffeq
                        twc = ex_to<lst>(twcex.subs((*imagepartSolu.begin()).op(0).lhs() == (*imagepartSolu.begin()).op(0).rhs()));
                        phase = ex_to<lst>(phaseex.subs((*imagepartSolu.begin()).op(0).lhs() == (*imagepartSolu.begin()).op(0).rhs()));
                        tw_coordiPhase = tw_coordiPhase.subs((*imagepartSolu.begin()).op(0).lhs() == (*imagepartSolu.begin()).op(0).rhs());
                        tw_coordi = tw_coordi.subs((*imagepartSolu.begin()).op(0).lhs() == (*imagepartSolu.begin()).op(0).rhs());
                        temdiffeq = temdiffeq.subs((*imagepartSolu.begin()).op(0).lhs() == (*imagepartSolu.begin()).op(0).rhs());

                        lst temvariables;
                        for (auto it = variables.begin();it != variables.end();it++) // remove (*imagepartSolu.begin()).op(0).lhs() from variables
                        {
                            if(*it != (*imagepartSolu.begin()).op(0).lhs())
                                temvariables.append(*it);
                        }
                        variables = temvariables;
                        (*imagepartSolu.begin()).op(0).lhs() = (*imagepartSolu.begin()).op(0).rhs();
                    }
                    else
                    {
                        cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        writetofile(solutions, dpndt_vars.op(0));

#ifdef GiNaCDE_gui
                        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
                        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
#endif // GiNaCDE_gui

                        return -1;
                    }

                }
                else
                {
                    exset temSymb = symbols(realpart); // checking present of symbols in realpart
                    if(temSymb.empty())
                    {
                        cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        writetofile(solutions, dpndt_vars.op(0));

                        #ifdef GiNaCDE_gui
                        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
                        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                        #endif // GiNaCDE_gui

                        return -1;
                    }

                    cond = false;
                    temex = collect_common_factors(realpart);
                    ex temexClt=_ex1;

                    if (is_a<mul>(temex))
                    {
                        for (size_t i = 0; i < nops(temex);i++)
                        {
                            if (!temex.op(i).has(U) && is_a<add>(temex.op(i)))
                            {
                                cond = true;
                                temexClt = temexClt*temex.op(i);
                            }
                        }
                    }
                    else if(!temex.has(U))
                    {
                        cond = true;
                    }
                    temex = temexClt;

                    if(cond)
                    {
                        temSymb = symbols(temex);
                        vector<ex> varSlct; // selecting parameter to get condition from real part
                        if(temSymb.size()>1)
                        {
                            std::vector<std::pair<ex, int>> temSymbMapWtDeg;
                            for (auto it = temSymb.begin();it != temSymb.end();it++)
                            {
                                if(is_polynomial(temex,*it))
                                    temSymbMapWtDeg.push_back(make_pair(*it,degree(temex,*it)));
                                else
                                    varSlct.push_back(*it);
                            }
                            if(temSymbMapWtDeg.size()>1)
                            {
                                sort(temSymbMapWtDeg.begin(), temSymbMapWtDeg.end(), [=](std::pair<ex, int>& a, std::pair<ex, int>& b){
                                    return a.second < b.second;});
                                varSlct.clear();
                                varSlct.push_back(temSymbMapWtDeg.begin()->first);
                            }
                        }
                        else if(!temSymb.empty())
                            varSlct.push_back(*temSymb.begin());

                        exsetlst realpartSolu;

                        if(!varSlct.empty())
                            realpartSolu = solve({temex},{*varSlct.begin()});
                        if(!realpartSolu.empty())
                        {
                            temdiffeq = imagepart;
                            remainingDiffpart = realpart;
                            solutions << "We derive solutions of imaginary part of Diff. Equ. with the condition:"
                                         "\n "<<(*realpartSolu.begin()).op(0).lhs()<<" == "<<(*realpartSolu.begin()).op(0).rhs()<<";"<<endl;

                            constraints.append((*realpartSolu.begin()).op(0).lhs()==(*realpartSolu.begin()).op(0).rhs());

                            ex twcex = twc,  phaseex = phase; // substitute (*realpartSolu.begin()).op(0).lhs() = (*realpartSolu.begin()).op(0).rhs() in twc, phase, temdiffeq
                            twc = ex_to<lst>(twcex.subs({(*realpartSolu.begin()).op(0).lhs() == (*realpartSolu.begin()).op(0).rhs()}));
                            phase = ex_to<lst>(phaseex.subs({(*realpartSolu.begin()).op(0).lhs() == (*realpartSolu.begin()).op(0).rhs()}));
                            tw_coordiPhase = tw_coordiPhase.subs({(*realpartSolu.begin()).op(0).lhs() == (*realpartSolu.begin()).op(0).rhs()});
                            tw_coordi = tw_coordi.subs({(*realpartSolu.begin()).op(0).lhs() == (*realpartSolu.begin()).op(0).rhs()});
                            temdiffeq = temdiffeq.subs({(*realpartSolu.begin()).op(0).lhs() == (*realpartSolu.begin()).op(0).rhs()});

                            lst temvariables;
                            for (auto it = variables.begin();it != variables.end();it++)// remove (*realpartSolu.begin()).op(0).lhs() from variables
                            {
                                if(*it != (*realpartSolu.begin()).op(0).lhs())
                                    temvariables.append(*it);
                            }
                            variables = temvariables;

                        }
                        else
                        {
                            cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                            solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                            writetofile(solutions, dpndt_vars.op(0));

#ifdef GiNaCDE_gui
                            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
                            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
#endif // GiNaCDE_gui

                            return -1;
                    }
                    }

                }
            }


            if(!cond && imagepart !=_ex0 && realpart != _ex0)
            {
                exmap coeffWtDpndtTermsR,coeffWtDpndtTermsI;
                collectAllCoeff(realpart,{U},true,coeffWtDpndtTermsR);
                collectAllCoeff(imagepart,{U},true,coeffWtDpndtTermsI);

                if(!coeffWtDpndtTermsR.empty()&&!coeffWtDpndtTermsI.empty())
                {
                    ex temRealPartActual=_ex0,temRealPartWtwild=_ex0, temImagePart;
                    exmap repls;

                    for(auto itr=coeffWtDpndtTermsR.begin(); itr!=coeffWtDpndtTermsR.end();itr++) // getting real part.
                    {
                        temRealPartActual=temRealPartActual+(itr->second)*(itr->first);
                    }


                    unsigned i = 0;
                    for(auto itr=coeffWtDpndtTermsR.begin(); itr!=coeffWtDpndtTermsR.end();itr++)
                    {
                        temRealPartWtwild=temRealPartWtwild+wild(i)*(itr->first);
                        i++;
                    }
                    for(auto itr=coeffWtDpndtTermsI.begin(); itr!=coeffWtDpndtTermsI.end();itr++) // getting imaginary part in rearranged form to match with imaginary part.
                    {
                        temImagePart=temImagePart+(itr->second)*(itr->first);
                    }


                    temImagePart.match(temRealPartWtwild,repls); // matching real part with imaginary part.

                    if(!repls.empty()) // comparison is successful
                    {
                        solutions << "Imaginary part and real part of Diff. Equ. becomes identical with following conditions: " << endl;
                        i=0;

                        for(auto itr=coeffWtDpndtTermsR.begin(); itr!=coeffWtDpndtTermsR.end();itr++)
                        {
                            solutions<<"("<<repls[wild(i)]<<")"<<"/"<<"("<<itr->second<<")";
                            constraints.append((repls[wild(i)])/(itr->second));
                            i++;

                            if(i<repls.size())
                                solutions<<" == ";
                            else if(i==repls.size())
                                solutions<<"; "<<endl;;
                        }
                        for(size_t constri=0; constri<nops(constraints)-1; constri++)
                        {
                            constraints[constri]=constraints[constri]==constraints[constri+1];
                        }
                        constraints.remove_last();

                        temdiffeq = temImagePart;
                        remainingDiffpart = temRealPartActual;
                        integrateDone=true;

                        solutions << "We derive solutions of imaginary part of Diff. Equ.." << endl;
                    }
                    else
                    {
                        cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
                        writetofile(solutions, dpndt_vars.op(0));

#ifdef GiNaCDE_gui
                        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
                        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
#endif // GiNaCDE_gui

                        return -1;
                    }


                }
                else
                {
                    temdiffeq=realpart+I*imagepart;
                    solutions << "The Diff. Equ. becomes: " << diffformchange(temdiffeq, (lst){U}, indpndt_vars) <<" = 0;" << endl;

                }

            }
            else if(!cond)
            {
                temdiffeq=realpart+I*imagepart;
                solutions << "The Diff. Equ. becomes: " << diffformchange(temdiffeq, (lst){U}, indpndt_vars) <<" = 0;" << endl;

            }
        }
        else
        {
            cout << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
            solutions << "Evaluation stop: unable to separate real and imaginary parts;" << endl;
            writetofile(solutions, dpndt_vars.op(0));

            #ifdef GiNaCDE_gui
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: Unable to separate real and imaginary parts;");
            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
            #endif // GiNaCDE_gui

            return -1;
        }
    }

/////////////////////////////////////////////////////////////////////////
    else
    {
        /// Checking integrability of ode
        ex odeprev;
        unsigned symnum = 0;
        bool cond = true;
        do
        {
            odeprev = integrate(temdiffeq, *indpndt_vars.begin());

            if(!odeprev.has(Integrate(wild(0), wild(1), wild(2))))
            {
                string tem1;
                ex tem2;

                if(symnum == 0)
                {
                    cout << "The Diff. Equ. is integrable;" << endl;
                    solutions << "The Diff. Equ. is integrable;" << endl;
                }
                symname = "ic_" + to_string(symnum+1); // number at subscript is number of integration.
                tem2 = reader(symname);

                #ifndef GiNaCDE_gui
                if(test)
                {
                    tem2 = 0;
                    solutions << "The assigned value(s) to the integration constant(s)-> " << endl;
                    solutions << "ic_1: " << 0 << endl;

                }
                else
                {
                    do
                    {
                        cout << "Do you assign a value to integration constant (ic_" <<symnum+1<< ")"<<"? ";
                        cin >> tem1;
                        if(tem1 == "h")
                        {
                            cout << "Type y for yes. \n"
                                    "Type n for no." ;
                        }
                        //cout << endl;
                    }while(tem1 != "y" && tem1 != "n");

                    if(tem1 == "y")
                    {
                       cout << "ic_" <<symnum+1<<": ";
                       cin >> tem1;
                       tem2 = reader(tem1);

                       if(symnum == 0)
                            solutions << "The assigned value(s) to the integration constant(s)-> " << endl;
                       solutions << "ic_" <<symnum+1<<": " << tem1 << endl;
                    }

                    if(!is_a<numeric>(tem2))
                        paraInDiffSolve.append(tem2);
                }

                #endif // GiNaCDE_gui

                temdiffeq = odeprev + tem2;

                symnum = symnum + 1;
            }
            else
                cond = false;

        }while(cond);

        temdiffeq = Simplify(expand(evaluate( temdiffeq )));

        #ifdef GiNaCDE_gui
        if(symnum)
        {
            string tem;

            string nlodeEq ;
            nlodeEq = "The Diff. Equ. is integrable;\nThe integration constant(s) are: ";

            GtkWidget *Eqlbl = gtk_label_new(&nlodeEq[0]);

            GtkWidget *table2 = gtk_grid_new();
            gtk_grid_attach(GTK_GRID(table2),Eqlbl,0,0,2,1);
            entrylbl.clear();
            entry.clear();
            entryText.clear();
            for(unsigned i=0; i<symnum; i++)
            {
                tem = "ic_"+to_string(i+1)+":";
                entrylbl.push_back(gtk_label_new(&tem[0]));
                entry.push_back(gtk_entry_new());
                gtk_grid_attach(GTK_GRID(table2),entrylbl[i],0,1+i,1,1);
                gtk_grid_attach(GTK_GRID(table2),entry[i],1,1+i,1,1);
                tem = "ic_"+to_string(i+1);
                gtk_entry_set_text(GTK_ENTRY(entry[i]), &tem[0]);
            }


            GtkWidget *dialog2=gtk_dialog_new();
            gtk_window_set_decorated(GTK_WINDOW(dialog2),FALSE);
            gtk_window_set_title(GTK_WINDOW(dialog2),"The integration constant(s) are:");
            gtk_window_set_modal( GTK_WINDOW(dialog2), TRUE );
            gtk_window_set_transient_for( GTK_WINDOW(dialog2), GTK_WINDOW(window));

            GtkAccelGroup *accel_group2=gtk_accel_group_new();
            GtkWidget *okbutton2=gtk_button_new_with_label("Ok");
            GtkWidget *cancelbutton2=gtk_button_new_with_label("Cancel");
            gtk_window_add_accel_group(GTK_WINDOW(dialog2),accel_group2);
            gtk_widget_add_accelerator(okbutton2,"clicked",accel_group2,GDK_KEY_Return,(GdkModifierType)0, (GtkAccelFlags)0);
            GtkWidget *hbox2=gtk_box_new(GTK_ORIENTATION_HORIZONTAL,10);
            gtk_box_pack_start(GTK_BOX(hbox2),okbutton2,FALSE,FALSE,0);
            gtk_box_pack_start(GTK_BOX(hbox2),cancelbutton2,FALSE,FALSE,0);
            gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog2)))),table2);
            //gtk_dialog_add_action_widget(GTK_DIALOG(dialog2), hbox2, GTK_RESPONSE_NONE);
            gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog2)))),hbox2);



            g_signal_connect(okbutton2,"clicked",G_CALLBACK(twc_okbutton_clicked),(gpointer)dialog2);
            g_signal_connect(cancelbutton2,"clicked",G_CALLBACK(cancelbutton_clicked),(gpointer)dialog2);
            for(unsigned i=0; i<entry.size(); i++)
                g_signal_connect(entry[i], "key-press-event", G_CALLBACK(key_press_cb),(gpointer)entrylbl[i]);
            gtk_widget_show_all(dialog2);
            gtk_dialog_run(GTK_DIALOG(dialog2));

          if(cancel)
          {
              cancel = false;

              #ifdef GiNaCDE_gui
              gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
              #endif // GiNaCDE_gui

              return 0;

          }

          unsigned onetime = 0;
          for(unsigned i = 0; i < symnum; i++)
          {
               ex tem1;
               tem1 = entryText[i];
               if(onetime == 0)
                    {solutions << "The assigned value(s) to the integration constant(s)-> " << endl;onetime = 1;}
               solutions << "ic_"<<i+1<<": " << tem1 << endl;;
               tem = "ic_"+to_string(i+1);
               ex icname = reader(tem);
               temdiffeq = Simplify(expand(temdiffeq.subs({icname==tem1},subs_options::algebraic)));
               if(!is_a<numeric>(tem1))
                    paraInDiffSolve.append(tem1);
          }
        }
        #endif // GiNaCDE_gui

        if((temdiffeq.subs(Diff(wild(0), wild(1), wild(2)) == symb_)).has(*indpndt_vars.begin())
            && temdiffeq.has(Diff(wild(0),wild(1),wild(2))))
        {
            cout << "Evaluation stop: independent variable present in diff. Equ.;" << endl;
            solutions << "Evaluation stop: independent variable present in diff. Equ.;" << endl;
            writetofile(solutions, dpndt_vars.op(0));

            #ifdef GiNaCDE_gui
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Evaluation stop: independent variable present in diff. Equ.;");
            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
            #endif // GiNaCDE_gui
            return -1;
        }

        if((tw_coordi != _ex0 || symnum != 0)&&!integrateDone)
        {
            solutions << "The Diff. Equ. becomes: " << diffformchange(temdiffeq, (lst){U}, indpndt_vars) <<" = 0;" << endl;
        }

    }


    if(tw_coordi != _ex0 && !phasepart)
    {
        if(output == mathematica)
            solutions << dpndt_vars.op(0) << " = " << U << "[xi], " << endl;
        else if(output == maple)
            solutions << dpndt_vars.op(0) << " = " << U << "(xi), " << endl;
        else
            solutions << dpndt_vars.op(0) << " = " << U << ", " << endl;
        solutions << "where xi = " << tw_coordi << ";" << endl;
    }
    else if( tw_coordi != _ex0 && phasepart )
    {
        if(output == mathematica)
            solutions << dpndt_vars.op(0) << " = " << U << "[xi]*exp(I*(" << tw_coordiPhase << "))," << endl;
        else if(output == maple)
            solutions << dpndt_vars.op(0) << " = " << U << "(xi)*exp(I*(" << tw_coordiPhase << "))," << endl;
        else
            solutions << dpndt_vars.op(0) << " = " << U << "*exp(I*(" << tw_coordiPhase << "))," << endl;
        solutions << "where U is the function of xi and xi = " << tw_coordi << ";" << endl;
    }

    if(!temdiffeq.has(Diff(wild(0),wild(1),wild(2))))
    {
        beginTime = chrono::high_resolution_clock::now();
        set<lst, ex_is_less>  solu_set_clt;

        lst coeffs_clt;

        solu_set_clt = solve({temdiffeq}, {U} );
        if(solu_set_clt.empty())
        {
            solutions<<"No solution exist;"<<endl;
            writetofile(solutions, dpndt_vars.op(0));

            #ifdef GiNaCDE_gui
            gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "No solution exist;");
            //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
            #endif // GiNaCDE_gui

            return -1;
        }

        if(solu_set_clt.size()==1)
            solutions <<"\n\n" << "The solution of the diff. equ. is:\n" << endl;
        else
            solutions <<"\n\n" << "The solutions of the diff. equ. are:\n" << endl;

        for(set<lst, ex_is_less>::const_iterator it = solu_set_clt.begin(); it != solu_set_clt.end(); it++)
        {

            if(tw_coordi != _ex0)
                solutions << dpndt_vars.op(0) << " = "  << Simplify(expand(((*it).op(0).rhs()*exp( I*tw_coordiPhase ))).subs(*indpndt_vars.begin()==tw_coordi)) <<";" << endl;
            else
                solutions << dpndt_vars.op(0) << " = "  << Simplify(expand((*it).op(0).rhs()*exp( I*tw_coordiPhase ))) <<";" << endl;

            solutions <<  "\n" << endl;
        }
        const std::chrono::time_point<std::chrono::system_clock> endTime = chrono::high_resolution_clock::now();
        auto dur = endTime-beginTime;
        cout << " Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;
        solutions << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" << endl;
        writetofile(solutions, dpndt_vars.op(0));

        #ifdef GiNaCDE_gui
        stringstream temstr;
        temstr << "The results are written in " << filename <<" file.             " << "Time:"
        << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()/1000.0 << " seconds" ;

        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, &temstr.str()[0]);
        //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
        #endif // GiNaCDE_gui

        return 0;
    }

    for(unsigned i = 0; i <paraInDiffSolve.nops(); i++)
    {
        variables.append(paraInDiffSolve.op(i));
    }
    paraInDiffSolve.remove_all();

    int ret = -1;
    if(method == F_expansion || method == mF_expansion)
    {
        ex Nvalue;

        if(NValue != 0)
        {
            if(NValue < 0)
            {
                solutions<<"Negative value of N is not allowed;"<<endl;
                writetofile(solutions, dpndt_vars.op(0));

                #ifdef GiNaCDE_gui
                gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Negative value of N is not allowed;");
                //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                #endif // GiNaCDE_gui

                return -1;
            }
            Nvalue = NValue;
            NValue = 0;
        }
        else
        {
            find_n_power find_n_power(0);
            find_n_power.n_pow(numer(temdiffeq), U);

            // Checking Nonlinearity of ODE //
            bool NonLinear = false;

            auto n_pow_clt_it = find_n_power.n_pow_clt.begin();
            while(!NonLinear && n_pow_clt_it != find_n_power.n_pow_clt.end())
            {
                if(coeff(*n_pow_clt_it, _n) > 1)
                    NonLinear = true;
                else
                    n_pow_clt_it++;
            }

            if(NonLinear == false)
                Nvalue = 1;
            else
            {
                /* identifying the turning point in find_n_power.n_pow_clt */
                std::vector<double> n_pow_value;
                std::vector<double> n_pow_max;
                bool infletious = false;
                double difference = 0;
                int NvalueInt;

                //dorat dorat;
                ex denoma =2;
                ex differStore;
                do
                {
                    NvalueInt = 0;

                    do
                    {
                        for(auto it = find_n_power.n_pow_clt.begin(); it != find_n_power.n_pow_clt.end(); it++)
                        {
                            n_pow_value.push_back(ex_to<numeric>(evalf((*it).subs(_n == NvalueInt/denoma))).to_double());
                        }
                        n_pow_max.push_back(*max_element(n_pow_value.begin(), n_pow_value.end()));

                        if(NvalueInt == 1)
                            difference = n_pow_max[1] - n_pow_max[0];

                        if(NvalueInt > 1)
                            differStore = (difference - (n_pow_max[NvalueInt] - n_pow_max[NvalueInt-1]));

                        if(NvalueInt > 1 && ((differStore.info(info_flags::positive) && differStore > GiNaC::pow(10,-2)) ||
                                             (differStore.info(info_flags::negative) && -differStore > GiNaC::pow(10,-2))))
                            infletious = true;

                        NvalueInt = NvalueInt + 1;

                    }while(!infletious && NvalueInt < 14*denoma);

                    n_pow_value.clear();
                    n_pow_max.clear();
                        denoma++;

                }while(!infletious && denoma < 5);

                if(evalf((NvalueInt-_ex2)/(denoma-_ex1)) > 10.0)
                {
                    solutions << "Evaluation stop: The value of N exceeds maximum value (10).\n"
                                  "Please try other value of N assigning NValue. Ex: NValue = 2." <<endl;
                    writetofile(solutions, dpndt_vars.op(0));
                    cout << "Evaluation stop: The value of N exceeds maximum value (10).\n"
                             "Please try other value of N assigning NValue. Ex: NValue = 2." <<endl;

                    #ifdef GiNaCDE_gui
                    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0,"Evaluation stop: The value of N exceeds maximum value (10).\n"
                             "Please try other value of N assigning NValue. Ex: NValue = 2.");
                    //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                    #endif // GiNaCDE_gui

                    return -1;
                }

                if(evalf((NvalueInt-_ex2)/(denoma-_ex1)) < 0)
                {
                    solutions<<"Negative value of N is not allowed;"<<endl;
                    writetofile(solutions, dpndt_vars.op(0));

                    #ifdef GiNaCDE_gui
                    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Negative value of N is not allowed;");
                    //gtk_show_uri(gdk_screen_get_default(),&CurrentPath[0],GDK_CURRENT_TIME,NULL);
                    #endif // GiNaCDE_gui

                    return -1;
                }

                Nvalue = (NvalueInt-_ex2)/(denoma-_ex1);

            }

        }

        solutions << "The value of N is: " << Nvalue << ";" << endl;
        cout << "The value of N is: " << Nvalue << ";" << endl;

        beginTime = chrono::high_resolution_clock::now();
        F_expans F_expans;

        ret = F_expans(numer(temdiffeq), U, dpndt_vars.op(0), *indpndt_vars.begin(), twc, tw_coordi, phase, tw_coordiPhase, variables, solutions, Nvalue, method, diffDenomSolu, remainingDiffpart);

    }
    else if(method == FIM)
    {
        beginTime = chrono::high_resolution_clock::now();
        fim fim;
        ret = fim(numer(temdiffeq), U, dpndt_vars.op(0), *indpndt_vars.begin(), twc, tw_coordi, phase, tw_coordiPhase, variables, solutions, diffDenomSolu, remainingDiffpart);
    }

    depend.clear(U, xi);

#ifdef GiNaCDE_gui
    resultsinDialog(solutions, dpndt_vars.op(0));
#endif //GiNaCDE_gui

    return ret;
}


ex checkSolu(const string& diff_equ, const string& solutions, const string& algebraic_solutions, const string& solutions_conditions)
{
    string solutionsR = replacestring(solutions, "==","=");


    ex diff_equex = reader(diff_equ);
    lst solutionslst={}, algebraic_solutionslst={},solutions_conditionslst={};
    vector<string> temSplit;

    temSplit = split(solutionsR,'=');
    solutionslst.append(reader(temSplit[0])==reader(temSplit[1]));

    // it handles complex NLPDE
    //if(solutions_conditions != "")
    //{
        string diff_equStr = diff_equ;
        if(diff_equStr.find("conjugate"))
        {
            temSplit[0].erase(std::remove(temSplit[0].begin(), temSplit[0].end(), ' '), temSplit[0].end());
            diff_equex = reader(replacestring(diff_equ, "conjugate("+temSplit[0]+")",temSplit[0]+"c_"));
            replaceI replaceI;
            ex replaceIex = subs(replaceI(reader(temSplit[1])),symb_==-I);
            solutionslst.append(reader(temSplit[0]+"c_")==replaceIex);
        }
    //}

    if(algebraic_solutions != "")
    {
        string algebraic_solutionsR = replacestring(algebraic_solutions, "==","=");

        algebraic_solutionsR.erase(std::remove(algebraic_solutionsR.begin(),algebraic_solutionsR.end(),'{'),algebraic_solutionsR.end());
        algebraic_solutionsR.erase(std::remove(algebraic_solutionsR.begin(),algebraic_solutionsR.end(),'}'),algebraic_solutionsR.end());

        vector<string> algebraic_solutionsSplit = split(algebraic_solutionsR,',');
        if(sizeof (algebraic_solutionsSplit) != 0)
        {
            for(auto itr = algebraic_solutionsSplit.begin(); itr != algebraic_solutionsSplit.end(); itr++)
            {
                temSplit = split(*itr,'=');
                algebraic_solutionslst.append(reader(temSplit[0])==reader(temSplit[1]));
            }
        }
    }

    if(solutions_conditions != "")
    {
        string solutions_conditionsR = replacestring(solutions_conditions, "==","=");

        solutions_conditionsR.erase(std::remove(solutions_conditionsR.begin(),solutions_conditionsR.end(),'{'),solutions_conditionsR.end());
        solutions_conditionsR.erase(std::remove(solutions_conditionsR.begin(),solutions_conditionsR.end(),'}'),solutions_conditionsR.end());

        vector<string> solutions_conditionsSplit = split(solutions_conditionsR,',');
        if(sizeof (solutions_conditionsSplit) != 0)
        {
            for(auto itr = solutions_conditionsSplit.begin(); itr != solutions_conditionsSplit.end(); itr++)
            {
                temSplit = split(*itr,'=');
                solutions_conditionslst.append(reader(temSplit[0])==reader(temSplit[1]));
            }
        }
    }


    if(nops(algebraic_solutionslst) != 0)
        diff_equex = subs(diff_equex,algebraic_solutionslst);
    if(nops(solutions_conditionslst) != 0)
        diff_equex = subs(diff_equex,solutions_conditionslst);

    //cout<<evaluate(subs(diff_equex,solutionslst))<<endl;
    diff_equex = (fullsimplify(evaluate(subs(diff_equex,solutionslst)),FuncSimp));

    return diff_equex;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
