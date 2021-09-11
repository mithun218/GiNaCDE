
/** @file GiNaCDE_gui.cpp
 *
 *  A gtk+3 GUI tool for solving differential equations using GiNaCDE library.*/


#include "GiNaCDE.h"
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <glib.h>


using namespace std;
using namespace GiNaC;

GtkWidget *window, *vbox, *asolve,*positivePartWid,*negativePartWid,*NValueWid, *fileNameWid, *dpndtVar, *indpndtVar, *diffEq,
          *extraVar, *outputFormat, *solvingMethod, *status_bar, *NValueWidLbl, *dpndtVarLbl, *indpndtVarLbl, *diffEqLbl,
          *positivePartLbl,*negativePartLbl,*asolveLbl, *solvingMethodLbl, *extraVarLbl, *textview1, *textview2, *entryw;

//GtkItemFactory *item_factory;
//GdkColor gtk_color;



GtkWidget *nlodeDeg;
vector<GtkWidget*> Aentry, Albl;
vector<ex> AentryText;

bool Cancel = false;

string odedegs;
const gchar *strLbl;

static void solve_equations_clicked(GtkWidget* button,gpointer dialog)
{
    GtkTextBuffer *buffer;
    GtkTextIter start, last;

    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview2));
    gtk_text_buffer_get_bounds(buffer, &start,&last);
    gtk_text_buffer_delete(buffer, &start, &last);

    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview1));
    gtk_text_buffer_get_start_iter(buffer, &start);
    gtk_text_buffer_get_end_iter(buffer, &last);
    string sysequText = gtk_text_buffer_get_text(buffer, &start, &last, FALSE);
    string varText =gtk_entry_get_text(GTK_ENTRY(entryw));

    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview2));
    gtk_text_buffer_get_start_iter(buffer, &start);
    if(sysequText=="" || varText=="")
    {
        gtk_text_buffer_insert(buffer, &start, "Error",-1);
        return;
    }

    try
    {
        set<lst, ex_is_less>  solu_set_clt;
        sysequText = "{"+sysequText+"}";
        varText = "{"+varText+"}";
        solu_set_clt = solve(ex_to<lst>(reader(sysequText)), ex_to<lst>(reader(varText)));

        if(!solu_set_clt.empty())
        {
            ostringstream str;
            for(set<lst, ex_is_less>::const_iterator it = solu_set_clt.begin(); it != solu_set_clt.end(); it++)
            {
                str << *it << endl;
            }
            gtk_text_buffer_insert(buffer, &start, &str.str()[0],-1);
        }
        else
        {
            gtk_text_buffer_insert(buffer, &start, "{}",-1);
        }


    }
    catch(GiNaC::parse_error)
    {
        gtk_text_buffer_insert(buffer, &start, "Error",-1);
        return;
    }


}

static void entry_insert_integer(GtkEditable *entry, char *New, gint len, gpointer position, gpointer data)
 {
   gboolean  faulty = FALSE;
   for (int i = 0; i < len; i++)
   {
     if(!isdigit (New[i]))
       faulty = TRUE;
   }
   if (faulty)
   {
     g_signal_stop_emission_by_name((gpointer)entry, "insert-text");
   }


 }

static void deg_okbutton_clicked(GtkWidget* okbutton,gpointer dialog)
{
    odedegs = gtk_entry_get_text(GTK_ENTRY(nlodeDeg));
    
    if(odedegs!="")
        gtk_widget_destroy(GTK_WIDGET(dialog));
}

static void Acoeff_okbutton_clicked(GtkWidget* okbutton,gpointer dialog)
{
    bool isValid = true;
    unsigned i = 0;

    AentryText.clear();

    string odedegsStr = gtk_entry_get_text(GTK_ENTRY(Aentry[stoi(odedegs)]));

    if(odedegsStr=="")
    {
        isValid = false; 

        /*char *str = g_strdup_printf ("<span font=\"14\" color=\"red\">"
                               "<b>\t\tRed: %d</b>"
                             "</span>",
                             value);*/
        strLbl = g_strdup_printf ("<span color=\"red\">"
                               "%s"
                             "</span>",
                             gtk_label_get_text(GTK_LABEL(Albl[stoi(odedegs)])));
        gtk_label_set_markup(GTK_LABEL(Albl[stoi(odedegs)]), strLbl);
        //gtk_widget_override_color(Albl[stoi(odedegs)], GTK_STATE_FLAG_NORMAL, &gtk_color);
    }
    else if (reader(odedegsStr)==_ex0)
    {
        /*char *str = g_strdup_printf ("<span font=\"14\" color=\"red\">"
                       "<b>\t\tRed: %d</b>"
                     "</span>",
                     value);*/
        strLbl = g_strdup_printf ("<span color=\"red\">"
                               "%s"
                             "</span>",
                             gtk_label_get_text(GTK_LABEL(Albl[stoi(odedegs)])));
        gtk_label_set_markup(GTK_LABEL(Albl[stoi(odedegs)]), strLbl);
        //gtk_widget_override_color(Albl[stoi(odedegs)], GTK_STATE_FLAG_NORMAL, &gtk_color);
    }

    else
        {
        for(auto itr = Aentry.begin(); itr != Aentry.end(); itr++ )
        {
            try
            {
                AentryText.push_back(reader(gtk_entry_get_text(GTK_ENTRY(*itr))));
                i = i+1;
            }
            catch(GiNaC::parse_error)
            {
                isValid = false;
                strLbl = g_strdup_printf ("<span color=\"red\">"
                               "%s"
                             "</span>",
                             gtk_label_get_text(GTK_LABEL(Albl[i])));
                gtk_label_set_markup(GTK_LABEL(Albl[i]), strLbl);


                //gtk_widget_modify_fg(Albl[i], GTK_STATE_NORMAL, &gtk_color);
                i = i+1;

            }
        }
    }

    if(isValid)
    {
        gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}

static void cancelbutton_clicked(GtkWidget* cancelbutton,gpointer dialog)
{
    Cancel = true;
    gtk_widget_destroy(GTK_WIDGET(dialog));
    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
}

/*preparation of menu item array */
/*................................................*/

void solve_equations()
{
    GtkTextBuffer *buffer;
    GtkTextIter start, last;

    GtkWidget *dialog=gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(dialog),"Solve the system of nonlinear polynomial equations.");
    gtk_window_set_modal( GTK_WINDOW(dialog), TRUE );
    gtk_window_set_transient_for( GTK_WINDOW(dialog), GTK_WINDOW(window));
    gtk_widget_set_size_request(GTK_WIDGET(dialog),1000,600);

    textview1 = gtk_text_view_new();
    //PangoFontDescription *fontdesc = pango_font_description_from_string("monospace 12");
    //gtk_widget_modify_font(textview1,fontdesc);
    GtkWidget *swin1 = gtk_scrolled_window_new (NULL, NULL);
    gtk_container_set_border_width (GTK_CONTAINER (swin1), 5);
    gtk_container_add (GTK_CONTAINER (swin1), textview1);

    gtk_box_pack_start((GtkBox*)(GtkDialog*)(gtk_dialog_get_content_area(GTK_DIALOG(dialog))),swin1,TRUE,TRUE,0);
    gtk_widget_set_tooltip_text(textview1,"Please enter the system of equations here.");
    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview1));
    gtk_text_buffer_get_start_iter(buffer, &start);
    gtk_text_buffer_get_end_iter(buffer, &last);
    gtk_text_buffer_insert(buffer, &start, "x+y,\nx^2-2+z,\nz+x",-1);

    entryw = gtk_entry_new();
    //gtk_widget_modify_font(entryw,fontdesc);
    gtk_entry_set_text(GTK_ENTRY(entryw),"x,y,z");
    gtk_box_pack_start((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog))),entryw,TRUE,TRUE,1);
    gtk_widget_set_tooltip_text(entryw,"Please enter the solution variables here.");

    textview2 = gtk_text_view_new();
    //gtk_widget_modify_font(textview2,fontdesc);
    GtkWidget *swin2 = gtk_scrolled_window_new (NULL, NULL);
    gtk_container_set_border_width (GTK_CONTAINER (swin2), 5);
    gtk_container_add (GTK_CONTAINER (swin2), textview2);
    gtk_box_pack_start((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog))),swin2,TRUE,TRUE,0);
    gtk_widget_set_tooltip_text(textview2,"Solutions are shown here.");

    GtkWidget *evaluateBtn=gtk_button_new();
    //fontdesc = pango_font_description_from_string("monospace 20");
    GtkWidget *evaluateBtnlabel = gtk_label_new("Solve");
    //gtk_widget_modify_font(evaluateBtnlabel,fontdesc);
    gtk_container_add(GTK_CONTAINER(evaluateBtn), evaluateBtnlabel);
    gtk_box_pack_start((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog))),evaluateBtn,TRUE,TRUE,0);
    gtk_widget_set_tooltip_text(evaluateBtn,"Click here to solve the equations.");

    gtk_widget_show_all(dialog);
    g_signal_connect(evaluateBtn,"clicked",G_CALLBACK(solve_equations_clicked),(gpointer)dialog);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);


}

void about()
{
     GtkWidget *infodialog = gtk_message_dialog_new(GTK_WINDOW(window),
                                   GTK_DIALOG_DESTROY_WITH_PARENT,
                                   GTK_MESSAGE_OTHER,
                                   GTK_BUTTONS_OK,
                   "GiNaCDE GUI (V1.0.0) build with GiNaC 1.7.6,\n"
                   "A GUI tool for GiNaCDE library to solve NLPDE or NLODE using\n"
                   "F-expansion and First integral methods.\n"
                   "Distributed under the terms of the GNU GPLv3\n"
                   "(see <https://www.gnu.org/licenses/gpl-3.0>)\n\n"
                   "This software uses GTK+3 toolkit\n(https://www.gtk.org)");

    gtk_window_set_icon_name(GTK_WINDOW(infodialog), "custom_icon");

    gtk_dialog_run (GTK_DIALOG(infodialog));
    gtk_widget_destroy (infodialog);
}



static void on_changed_solvingMethod(GtkWidget *_solvingMethod, gpointer user_data)
{
    string activeText = gtk_combo_box_text_get_active_text (GTK_COMBO_BOX_TEXT(_solvingMethod));
    if(activeText != "FIM")
    {
        gtk_widget_set_sensitive(asolve,TRUE);
        gtk_widget_set_sensitive(positivePartWid,TRUE);
        gtk_widget_set_sensitive(negativePartWid,TRUE);
    }
    else
    {
        gtk_widget_set_sensitive(asolve,FALSE);
        gtk_widget_set_sensitive(positivePartWid,FALSE);
        gtk_widget_set_sensitive(negativePartWid,FALSE);
    }

    gtk_label_set_text(GTK_LABEL(solvingMethodLbl), gtk_label_get_text(GTK_LABEL(solvingMethodLbl)));
    //gtk_widget_modify_fg(solvingMethodLbl, GTK_STATE_NORMAL, NULL);

}


static gboolean key_press_cb(GtkEditable *w, gpointer data)
{
    gtk_label_set_text(GTK_LABEL(data), gtk_label_get_text(GTK_LABEL(data)));
    //gtk_widget_modify_fg(GTK_WIDGET(data), GTK_STATE_NORMAL, NULL);
    return FALSE;
}

static gboolean nlodeDeg_key_press_cb(GtkEditable *w, GdkEvent *ev, gpointer data)
{
    GdkEventKey *key = (GdkEventKey*)ev;
    if(key->keyval == GDK_KEY_Escape) // preventing closing of dialog for escape key
    {
        return TRUE;
    }
    return FALSE;
}

static gboolean Aentry_key_press_cb(GtkEditable *w, GdkEvent *ev, gpointer data)
{
    GdkEventKey *key = (GdkEventKey*)ev;
    if(key->keyval == GDK_KEY_Escape) // preventing closing of dialog for escape key
    {
        return TRUE;
    }
    gtk_label_set_text(GTK_LABEL(data), gtk_label_get_text(GTK_LABEL(data)));
    //gtk_widget_modify_fg(GTK_WIDGET(data), GTK_STATE_NORMAL, NULL);
    return FALSE;
}

static int on_clicked_evaluatebtn(GtkWidget *_evaluatebtn, gpointer user_data)
{
    gtk_statusbar_push(GTK_STATUSBAR(status_bar), 0, "Evaluating...");

    const string text_dpndtVar = gtk_entry_get_text(GTK_ENTRY(dpndtVar)),
                 text_indpndtVar = gtk_entry_get_text(GTK_ENTRY(indpndtVar)),
                 text_diffEq = gtk_entry_get_text(GTK_ENTRY(diffEq)),
                 text_solvingMethod = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(solvingMethod)),
                 text_extraVar = gtk_entry_get_text(GTK_ENTRY(extraVar)),
                 text_NValueWid = gtk_entry_get_text(GTK_ENTRY(NValueWid)),
                 text_fileNameWid = gtk_entry_get_text(GTK_ENTRY(fileNameWid));

    vector<string> temindVars;

    if(text_fileNameWid != "")
        filename = text_fileNameWid;

    bool isValid = true;
    ex  inputdiff,dpndt_var;
    lst indpndt_vars;


    if(text_solvingMethod == "")
    {

        strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(solvingMethodLbl)));
        gtk_label_set_markup(GTK_LABEL(solvingMethodLbl), strLbl);
        //gtk_widget_modify_fg(solvingMethodLbl, GTK_STATE_NORMAL, &gtk_color);
        isValid = false;
    }


    if(text_dpndtVar == "" || !is_a<symbol>(reader(text_dpndtVar)) )
    {
        strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(dpndtVarLbl)));
        gtk_label_set_markup(GTK_LABEL(dpndtVarLbl), strLbl);        
        //gtk_widget_modify_fg(dpndtVarLbl, GTK_STATE_NORMAL, &gtk_color);
        isValid = false;
    }
    else
    {
        try
        {
            dpndt_var = reader(text_dpndtVar);
            depend.clear(dpndt_var);
        }
        catch(GiNaC::parse_error)
        {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                       "%s"
                     "</span>",
                     gtk_label_get_text(GTK_LABEL(dpndtVarLbl)));
            gtk_label_set_markup(GTK_LABEL(dpndtVarLbl), strLbl);              
            //gtk_widget_modify_fg(dpndtVarLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
        }
    }

    if(text_indpndtVar == "")
    {
        strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(indpndtVarLbl)));
        gtk_label_set_markup(GTK_LABEL(indpndtVarLbl), strLbl);          
        //gtk_widget_modify_fg(indpndtVarLbl, GTK_STATE_NORMAL, &gtk_color);
        isValid = false;
    }
    else
    {
        try
        {
            if(text_indpndtVar.find(',')!=string::npos)
            {
                temindVars.clear();
                temindVars = split(text_indpndtVar,',');

                for(unsigned i=0; i<temindVars.size(); i++)
                    indpndt_vars.append(reader(temindVars[i]));
            }
            else
            {
                indpndt_vars.append(reader(text_indpndtVar));
            }
            for(auto itr = indpndt_vars.begin(); itr != indpndt_vars.end(); itr++)
            {
                if(!is_a<symbol>(*itr))
                {
                    strLbl = g_strdup_printf ("<span color=\"red\">"
                           "%s"
                         "</span>",
                         gtk_label_get_text(GTK_LABEL(indpndtVarLbl)));
                    gtk_label_set_markup(GTK_LABEL(indpndtVarLbl), strLbl);                 
                    //gtk_widget_modify_fg(indpndtVarLbl, GTK_STATE_NORMAL, &gtk_color);
                    isValid = false;
                }
            }
            if(isValid)
            {
                for(auto it = indpndt_vars.begin(); it != indpndt_vars.end(); it++)
                {
                    depend(dpndt_var, {*it});
                }
            }
        }
        catch(GiNaC::parse_error)
        {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                       "%s"
                     "</span>",
                     gtk_label_get_text(GTK_LABEL(indpndtVarLbl)));
            gtk_label_set_markup(GTK_LABEL(indpndtVarLbl), strLbl);             
            //gtk_widget_modify_fg(indpndtVarLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
        }
    }

    if(text_extraVar != "")
    {
        try
        {
            if(text_extraVar.find(',')!=string::npos)
            {
                temindVars.clear();
                temindVars = split(text_extraVar,',');

                for(unsigned i=0; i<temindVars.size(); i++)
                    paraInDiffSolve.append(reader(temindVars[i]));
            }
            else
            {
                paraInDiffSolve.append(reader(text_extraVar));
            }

            const exset symClt = symbols(paraInDiffSolve);
            for(auto itr = symClt.begin(); itr != symClt.end(); itr++)
            {
                if(!is_a<symbol>(*itr))
                {

                    strLbl = g_strdup_printf ("<span color=\"red\">"
                           "%s"
                         "</span>",
                         gtk_label_get_text(GTK_LABEL(extraVarLbl)));
                    gtk_label_set_markup(GTK_LABEL(extraVarLbl), strLbl);                 
                    //gtk_widget_modify_fg(extraVarLbl, GTK_STATE_NORMAL, &gtk_color);
                    isValid = false;
                }
            }

        }
        catch(GiNaC::parse_error)
        {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(extraVarLbl)));
            gtk_label_set_markup(GTK_LABEL(extraVarLbl), strLbl);             
            //gtk_widget_modify_fg(extraVarLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
        }
    }


    if(text_diffEq == "")
    {
        strLbl = g_strdup_printf ("<span color=\"red\">"
               "%s"
             "</span>",
             gtk_label_get_text(GTK_LABEL(diffEqLbl)));
        gtk_label_set_markup(GTK_LABEL(diffEqLbl), strLbl);         
        //gtk_widget_modify_fg(diffEqLbl, GTK_STATE_NORMAL, &gtk_color);
        isValid = false;
    }
    else
    {
        try
        {
            inputdiff = reader(text_diffEq);
        }
        catch(GiNaC::parse_error)
        {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(diffEqLbl)));
            gtk_label_set_markup(GTK_LABEL(diffEqLbl), strLbl);             
            //gtk_widget_modify_fg(diffEqLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
        }
    }

    if(text_NValueWid == "")
    {
        NValue = 0; // for auto-evaluation of N
    }
    else
    {
        try
        {
            NValue = reader(text_NValueWid);
        }
        catch(GiNaC::parse_error)
        {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(NValueWidLbl)));
            gtk_label_set_markup(GTK_LABEL(NValueWidLbl), strLbl);             
            //gtk_widget_modify_fg(NValueWidLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
        }

        if(!(NValue).info(info_flags::rational)) // numeric value is only allowed
          {
            strLbl = g_strdup_printf ("<span color=\"red\">"
                   "%s"
                 "</span>",
                 gtk_label_get_text(GTK_LABEL(NValueWidLbl)));
            gtk_label_set_markup(GTK_LABEL(NValueWidLbl), strLbl);             
            //gtk_widget_modify_fg(NValueWidLbl, GTK_STATE_NORMAL, &gtk_color);
            isValid = false;
          }
    }


    if(!isValid)
    {
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");
        return 0;
    }



    const string text_outputFormat = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(outputFormat)),
                 text_positivePart = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(positivePartWid)),
                 text_negativePart = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(negativePartWid)),
                 text_asolve = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(asolve));
    int method;

    if(text_positivePart == "yes")
        positivePart = true;
    else if(text_positivePart == "no")
        positivePart = false;


    if(text_negativePart == "yes")
        negativePart = true;
    else if(text_negativePart == "no")
        negativePart = false;

    if(!positivePart && !negativePart)
    {
        gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Minimum one part (positivePart or negativePart) should be \"yes\";");
        return 0;
    }



    if(text_solvingMethod == "FIM")
        method = FIM;
    else if(text_solvingMethod == "F_expansion")
    {
        method = F_expansion;
    }
    else if(text_solvingMethod == "mF_expansion")
    {
        method = mF_expansion;
    }

    if( method != FIM)
    {
        degAcoeff.remove_all();

        GtkWidget *nlodeDegLbl=gtk_label_new((gchar*)"Highest positive integer (delta)\n of 1st order NLODE (A.E.)");
        nlodeDeg=gtk_entry_new();

        GtkWidget *table = gtk_grid_new();
        gtk_grid_attach(GTK_GRID(table),nlodeDegLbl,0,0,1,1);
        gtk_grid_attach(GTK_GRID(table),nlodeDeg,1,0,1,1);

        GtkWidget *dialog=gtk_dialog_new();
        gtk_window_set_decorated(GTK_WINDOW(dialog),FALSE);
        gtk_window_set_title(GTK_WINDOW(dialog),"Enter delta:");
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
        gtk_container_add(GTK_CONTAINER ((GtkBox *) (GtkDialog *) (gtk_dialog_get_content_area(GTK_DIALOG(dialog)))),hbox);


        g_signal_connect(okbutton,"clicked",G_CALLBACK(deg_okbutton_clicked),(gpointer)dialog);
        g_signal_connect(nlodeDeg, "insert-text", G_CALLBACK(entry_insert_integer), NULL);
        g_signal_connect(cancelbutton,"clicked",G_CALLBACK(cancelbutton_clicked),(gpointer)dialog);
        g_signal_connect(nlodeDeg,"key-press-event",G_CALLBACK(nlodeDeg_key_press_cb),(gpointer)dialog); // preventing closing of dialog for escape key
        //gtk_window_set_icon_name(GTK_WINDOW(dialog), "custom_icon");
        gtk_widget_show_all(dialog);
        gtk_dialog_run(GTK_DIALOG(dialog));

        if(Cancel)
            {Cancel = false;gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready"); gtk_widget_set_sensitive(window,TRUE);return 0;}

        string tem;
        int int_nlodeDeg = stoi(odedegs);
        degAcoeff.append(reader(odedegs));

        string nlodeEq = "A_0";
        for(int i = 1; i <= int_nlodeDeg; i++)
        {
            nlodeEq = nlodeEq + "+A_" + to_string(i) + "*F^" + to_string(i);
        }

        if(method == mF_expansion)
            nlodeEq = "1st order NLODE (general):\n F' = " + nlodeEq;
        else
            nlodeEq = "1st order NLODE (general):\n F' = sqrt( "+ nlodeEq+" ) ";

        GtkWidget *Eqlbl = gtk_label_new(&nlodeEq[0]);
        table = gtk_grid_new();
        gtk_grid_attach(GTK_GRID(table),Eqlbl,0,0,2,1);
        Aentry.clear();
        Albl.clear();
        for(int i=0; i<=int_nlodeDeg; i++)
        {
            tem = "A_"+to_string(i)+":";
            Albl.push_back(gtk_label_new(&tem[0]));
            Aentry.push_back(gtk_entry_new());
            const string symname = "A_"+to_string(i);
            gtk_entry_set_text(GTK_ENTRY(Aentry[i]),&symname[0]);
            gtk_grid_attach(GTK_GRID(table),Albl[i],0,1+i,1,1);
            gtk_grid_attach(GTK_GRID(table),Aentry[i],1,1+i,1,1);
        }

        dialog=gtk_dialog_new();
        gtk_window_set_decorated(GTK_WINDOW(dialog),FALSE);
        gtk_window_set_title(GTK_WINDOW(dialog),"Enter coefficients:");
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


        g_signal_connect(okbutton,"clicked",G_CALLBACK(Acoeff_okbutton_clicked),(gpointer)dialog);
        g_signal_connect(cancelbutton,"clicked",G_CALLBACK(cancelbutton_clicked),(gpointer)dialog);
        for(unsigned i=0; i<Aentry.size(); i++)
            g_signal_connect(Aentry[i], "key-press-event", G_CALLBACK(Aentry_key_press_cb),(gpointer)Albl[i]); // preventing closing of dialog for escape key
        //gtk_window_set_icon_name(GTK_WINDOW(dialog), "custom_icon");
        gtk_widget_show_all(dialog);
        gtk_dialog_run(GTK_DIALOG(dialog));

        if(Cancel)
            {Cancel = false;gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");    gtk_widget_set_sensitive(window,TRUE);return 0;}

        for(int i = 0; i <= int_nlodeDeg; i++)
        {
           degAcoeff.append(AentryText[i]);
        }

    }

    g_free((gpointer)strLbl);

    if(text_asolve == "yes")
        ASolve = true;
    else if(text_asolve == "no")
        ASolve = false;



    if(text_outputFormat == "maple")
        output = maple;
    else
        output = mathematica;

    twcPhase={lst{},lst{}};
    //getting current directory path
    char cCurrentPath[FILENAME_MAX];
    GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
    CurrentPath = cCurrentPath;
    #ifdef WINDOWS
        CurrentPath = CurrentPath + "\\" + filename;
    #else
        CurrentPath = CurrentPath + "/" + filename;
    #endif

    desolve(inputdiff, {dpndt_var}, method);

    return 0;
}



int main(int argc,char* argv[])
{

    gtk_init (&argc, &argv);

    GtkWidget *menu_bar, *outputFormatLbl, *fileNameLbl, *evaluateBtn, *table;

    //gdk_color_parse ("#FF0000", &gtk_color);
    PangoFontDescription *fontdesc = pango_font_description_from_string("monospace 10");

    status_bar = gtk_statusbar_new ();
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window),"GiNaCDE GUI!");
    g_signal_connect (G_OBJECT (window), "destroy",G_CALLBACK (gtk_main_quit), NULL);
    gtk_container_set_border_width (GTK_CONTAINER (window), 10);
    gtk_widget_set_size_request(GTK_WIDGET(window),800,500);
    gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);

    vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,1);
    gtk_container_add(GTK_CONTAINER(window),vbox);
    gtk_widget_show(vbox);



    menu_bar =gtk_menu_bar_new ();

    GtkWidget* menuitem1 = gtk_menu_item_new_with_mnemonic ("_File");
    GtkWidget* menuitem1_1 = gtk_menu_new ();
    GtkWidget* item_Quit= gtk_menu_item_new_with_label ("Quit");
    gtk_menu_shell_append (GTK_MENU_SHELL (menuitem1_1), item_Quit);
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (menuitem1), menuitem1_1);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu_bar), menuitem1);

    GtkWidget* menuitem2 = gtk_menu_item_new_with_mnemonic ("_Solve");
    GtkWidget* menuitem2_1 = gtk_menu_new ();
    GtkWidget* item_Solve= gtk_menu_item_new_with_label ("Solve system of equations");
    gtk_menu_shell_append (GTK_MENU_SHELL (menuitem2_1), item_Solve);
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (menuitem2), menuitem2_1);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu_bar), menuitem2);

    GtkWidget* menuitem3 = gtk_menu_item_new_with_mnemonic ("_About");
    GtkWidget* menuitem3_1 = gtk_menu_new ();
    GtkWidget* item_GiNaC= gtk_menu_item_new_with_label ("GiNaCDE GUI");
    gtk_menu_shell_append (GTK_MENU_SHELL (menuitem3_1), item_GiNaC);
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (menuitem3), menuitem3_1);
    gtk_menu_shell_append (GTK_MENU_SHELL (menu_bar), menuitem3);

    gtk_box_pack_start(GTK_BOX(vbox),menu_bar,FALSE,TRUE,0);
    gtk_widget_show(menu_bar);



    dpndtVarLbl=gtk_label_new((gchar*)"Dependent variable");
    indpndtVarLbl=gtk_label_new((gchar*)"Independent variable(s)");
    diffEqLbl=gtk_label_new((gchar*)"Diff. Equation");
    outputFormatLbl=gtk_label_new((gchar*)"Output format");
    NValueWidLbl=gtk_label_new((gchar*)"Value of N");
    fileNameLbl=gtk_label_new((gchar*)"Name of file");
    solvingMethodLbl=gtk_label_new((gchar*)"Solving method");
    positivePartLbl=gtk_label_new((gchar*)"positivePart");
    negativePartLbl=gtk_label_new((gchar*)"negativePart");
    asolveLbl=gtk_label_new((gchar*)"ASolve");
    extraVarLbl=gtk_label_new((gchar*)"Supply extra constant parameters\npresent in diff. equ.");

    dpndtVar=gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(dpndtVar),(gchar*)"u");
    indpndtVar=gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(indpndtVar),(gchar*)"t,x");
    diffEq=gtk_entry_new();
    //gtk_widget_modify_font(diffEq,fontdesc);
    gtk_entry_set_text(GTK_ENTRY(diffEq),(gchar*)"p*Diff(u,t,2)+q*u");
    NValueWid=gtk_entry_new();
    fileNameWid=gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(fileNameWid),(gchar*)"outputResults.txt");
    extraVar=gtk_entry_new();

    outputFormat=gtk_combo_box_text_new ();
    solvingMethod=gtk_combo_box_text_new ();
    positivePartWid=gtk_combo_box_text_new ();
    negativePartWid=gtk_combo_box_text_new ();
    asolve=gtk_combo_box_text_new ();

    /*inserting text in combo box text*/
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(outputFormat), (gchar*)"maple");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(outputFormat), (gchar*)"mathematica");
    gtk_combo_box_set_active (GTK_COMBO_BOX(outputFormat), 0);

    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(solvingMethod), (gchar*)"FIM");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(solvingMethod), (gchar*)"F_expansion");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(solvingMethod), (gchar*)"mF_expansion");
    gtk_combo_box_set_active(GTK_COMBO_BOX(solvingMethod), 0);

    gtk_widget_set_sensitive(positivePartWid,FALSE);
    gtk_widget_set_sensitive(negativePartWid,FALSE);
    gtk_widget_set_sensitive(asolve,FALSE);

    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(positivePartWid), (gchar*)"yes");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(positivePartWid), (gchar*)"no");
    gtk_combo_box_set_active (GTK_COMBO_BOX(positivePartWid), 0);

    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(negativePartWid), (gchar*)"yes");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(negativePartWid), (gchar*)"no");
    gtk_combo_box_set_active (GTK_COMBO_BOX(negativePartWid), 0);

    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(asolve), (gchar*)"no");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(asolve), (gchar*)"yes");
    gtk_combo_box_set_active (GTK_COMBO_BOX(asolve), 0);

    evaluateBtn=gtk_button_new_with_label("Evaluate");

    /* tooltip text */
    gtk_widget_set_tooltip_text(dpndtVar,"Please provide dependent variable of input diff equ.");
    gtk_widget_set_tooltip_text(indpndtVar,"Please provide independent variable(s) of input diff equ. "
                            "Example: t,x,y,z");
    gtk_widget_set_tooltip_text(diffEq,"Examples of input diffs:\n1. Diff(u,t,1)-u*u*a*Diff(u,x,1)+Diff(u,x,3), \n"
                                                                 "2. Diff(u,t,1) + Diff(u,x,1) + u*Diff(u,x,1) - Diff(Diff(u,x,2),t,1), \n"
                                                                 "3. I*Diff(u,t,1) + p*Diff(u,x,2) + q*u*u*conjugate(u)");
    gtk_widget_set_tooltip_text(outputFormat,"Output results are saved in maple or mathematica format.");
    gtk_widget_set_tooltip_text(NValueWid,"For FIM, 'N' represents number of terms in sum(a_i(X)*Y^i,i=0..N) where X=u(x), Y=diff(u(x),x). For FIM, N = 1 and N = 2 are only allowed. "
                            "For F-expansion and modified F-expansion methods, the solutions of input NLPDE is expressed by a finite power series where presents 'N+1' terms.\n"
                            "By default, for FIM, N=1 and for F-expansion, modified F-expansion methods 'N' is evaluated automatically if N is not assigned to any value.");
    gtk_widget_set_tooltip_text(solvingMethod,"FIM for first integral method.\n"
                                              "F_expansion for F-expansion method.\n"
                                              "mF_expansion for modified F-expansion method.");
    gtk_widget_set_tooltip_text(asolve,"yes for solving algebraic system for parameters A_i (i=0,1,..delta) in A.E.., "
                                       "otherwise choose no.");
    gtk_widget_set_tooltip_text(extraVar,"GiNaCDE library solve the generated overdetermined system of algebraic equations only for "
                                         "internally generated constant parameters. If you wish to solve for any constant "
                                         "parameters appeared in differential equation provide them here.\n Example: p,q");

    gtk_widget_set_tooltip_text(fileNameWid,"Supply the file name where output results will be saved.");
    gtk_widget_set_tooltip_text(positivePartWid,"Choose yes to take only positive part in finite power series,"
                                                "otherwise choose no to take negative part only.");
    gtk_widget_set_tooltip_text(negativePartWid,"Choose yes to take only negative part in finite power series,"
                                                "otherwise choose no to take positive part only.");

    gtk_widget_set_tooltip_text(evaluateBtn,"Click here to start the evaluations.");

    status_bar = gtk_statusbar_new ();

    
    /*Grid creation and grid attach with label,entry,hbox widget*/
    table=gtk_grid_new();
    gtk_grid_set_row_spacing (GTK_GRID (table),15);
    gtk_grid_set_column_spacing(GTK_GRID(table),5);
    gtk_container_set_border_width (GTK_CONTAINER(table), 5);


    gtk_grid_attach(GTK_GRID(table),dpndtVarLbl,0,0,1,1);
    gtk_grid_attach(GTK_GRID(table),dpndtVar,1,0,1,1);
    gtk_grid_attach(GTK_GRID(table),indpndtVarLbl,2,0,1,1);
    gtk_grid_attach(GTK_GRID(table),indpndtVar,3,0,1,1);

    gtk_grid_attach(GTK_GRID(table),diffEqLbl,0,1,1,1);
    gtk_grid_attach(GTK_GRID(table),diffEq,1,1,3,1);

    gtk_grid_attach(GTK_GRID(table),outputFormatLbl,0,2,1,1);
    gtk_grid_attach(GTK_GRID(table),outputFormat,1,2,1,1);
    gtk_grid_attach(GTK_GRID(table),NValueWidLbl,2,2,1,1);
    gtk_grid_attach(GTK_GRID(table),NValueWid,3,2,1,1);

    gtk_grid_attach(GTK_GRID(table),fileNameLbl,0,3,1,1);
    gtk_grid_attach(GTK_GRID(table),fileNameWid,1,3,1,1);
    gtk_grid_attach(GTK_GRID(table),solvingMethodLbl,2,3,1,1);
    gtk_grid_attach(GTK_GRID(table),solvingMethod,3,3,1,1);

    gtk_grid_attach(GTK_GRID(table),positivePartLbl,0,4,1,1);
    gtk_grid_attach(GTK_GRID(table),positivePartWid,1,4,1,1);
    gtk_grid_attach(GTK_GRID(table),negativePartLbl,2,4,1,1);
    gtk_grid_attach(GTK_GRID(table),negativePartWid,3,4,1,1);

    gtk_grid_attach(GTK_GRID(table),asolveLbl,0,5,1,1);
    gtk_grid_attach(GTK_GRID(table),asolve,1,5,1,1);

    gtk_grid_attach(GTK_GRID(table),extraVarLbl,0,6,1,1);
    gtk_grid_attach(GTK_GRID(table),extraVar,1,6,1,1);

    gtk_grid_attach(GTK_GRID(table),evaluateBtn,3,7,1,1);

    gtk_grid_attach(GTK_GRID(table),status_bar,0,8,2,2);


    gtk_box_pack_start(GTK_BOX(vbox),table,FALSE,FALSE,0);

    gtk_statusbar_push (GTK_STATUSBAR(status_bar), 0, "Ready");


    /* connecting signals in menu items by signal callback */
    g_signal_connect_swapped (item_Quit, "activate", G_CALLBACK (gtk_main_quit), window);
    g_signal_connect (item_Solve, "activate", G_CALLBACK (solve_equations), NULL);
    g_signal_connect (item_GiNaC, "activate", G_CALLBACK (about), NULL);

    /* connecting signals in text entry field by signal callback */
    g_signal_connect(solvingMethod, "changed", G_CALLBACK (on_changed_solvingMethod), NULL);
    g_signal_connect(evaluateBtn,"clicked",G_CALLBACK(on_clicked_evaluatebtn),NULL);
    g_signal_connect(NValueWid, "changed", G_CALLBACK(key_press_cb), (gpointer)NValueWidLbl);
    g_signal_connect(dpndtVar, "changed", G_CALLBACK(key_press_cb),(gpointer)dpndtVarLbl);
    g_signal_connect(indpndtVar, "changed", G_CALLBACK(key_press_cb),(gpointer)indpndtVarLbl);
    g_signal_connect(diffEq, "changed", G_CALLBACK(key_press_cb),(gpointer)diffEqLbl);
    g_signal_connect(extraVar, "changed", G_CALLBACK(key_press_cb),(gpointer)extraVarLbl);

    gtk_widget_show_all(window);
    gtk_main();

     return 0;
  }
