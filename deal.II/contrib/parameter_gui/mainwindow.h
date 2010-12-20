
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTreeWidget>
#include <QDialog>
#include <QSettings>

#include "info_message.h"
#include "parameter_delegate.h"
#include "xml_parameter_reader.h"
#include "xml_parameter_writer.h"

/**
 * The MainWindow class of the the parameterGUI.
 * The parameterGUI is a graphical user interface for editing parameter files based on the XML format,
 * created by the ParameterHandler::print_parameters() function with ParameterHandler::XML as second argument.
 * Please see @ref Representation of Parameters in the documentation of the @ref ParameterHandler class for more details.
 * The MainWindow class provides the basic functionality of the GUI as save- and load-file-actions and documentation.
 * The parameterGUI provides special editors for the different types of parameters defined in the ParameterHandler class.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
 */
class MainWindow : public QMainWindow
{
  Q_OBJECT

  public:
				     /**
				      * Constructor.
				      * If a <tt>filename</tt> is given,
				      * the MainWindow tries to open
				      * and parse the file.
				      */
    MainWindow(const QString  &filename = "");

  protected:
				     /**
				      * Reimplemented from QMainWindow.
				      * We ask, if changes should be saved.
				      */
    void closeEvent(QCloseEvent *event);

  private slots:

				     /**
				      * Open a parameter file.
				      */
    void open();
				     /**
				      * Save the parameter file.
				      */
    bool save();
				     /**
				      * Open a file dialog to save the parameter file.
				      */
    bool save_as();
				     /**
				      * Show some information on the parameterGUI
				      */
    void about();

				     /**
				      * A <tt>slot</tt> that should be always called,
				      * if parameter values were changed.
				      */
    void tree_was_modified();

  private:
				     /**
				      * Show an information dialog, how
				      * parameters can be edited.
				      */
    void show_message ();
				     /**
				      * This function creates all actions.
				      */
    void create_actions();
				     /**
				      * This function creates all menus.
				      */
    void create_menus();
				     /**
				      * This function checks, if parameters were changed
				      * and show a dialog, if changes should be saved.
				      * This function should be always called,
				      * before open a new parameter file or before closing the GUI
				      */
    bool maybe_save ();
				     /**
				      * Save parameters to <tt>filename</tt> in XML format.
				      */
    bool save_file (const QString &filename);
				     /**
				      * Load parameters from <tt>filename</tt> in XML format.
				      */
    void load_file (const QString &filename);
				     /**
				      * This functions writes the current <tt>filename</tt> to the window title.
				      */
    void set_current_file (const QString  &filename);
				     /**
				      * This is the tree structure in which we store all parameters.
				      */
    QTreeWidget * tree_widget;
				     /**
				      * This menu provides all file actions as <tt>open</tt>, <tt>save</tt>, <tt>save as</tt>
				      * and <tt>exit</tt>
				      */
    QMenu * file_menu;
				     /**
				      * This menu provides some informations <tt>about</tt> the parameterGUI
				      * and <tt>about Qt</tt>
				      */
    QMenu * help_menu;
				     /**
				      * QAction <tt>open</tt> a file.
				      */
    QAction * open_act;
				     /**
				      * QAction <tt>save</tt> a file.
				      */
    QAction * save_act;
				     /**
				      * QAction <tt>save as</tt> a file.
				      */
    QAction * save_as_act;
				     /**
				      * QAction <tt>exit</tt> the GUI.
				      */
    QAction * exit_act;
				     /**
				      * QAction <tt>about</tt> the parameterGUI.
				      */
    QAction * about_act;
				     /**
				      * QAction <tt>about</tt> Qt.
				      */
    QAction * about_qt_act;
				     /**
				      * This value stores the current <tt>filename</tt> we work on.
				      */
    QString  current_file;
				     /**
				      * This dialog shows a short information message after loading a file.
				      */
    InfoMessage * info_message;
				     /**
				      * An object for storing user settings.
				      */
    QSettings * gui_settings;
};


#endif
