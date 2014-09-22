// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by Martin Steigemann and Wolfgang Bangerth
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <QtGui>

#include "mainwindow.h"
#include "parameter_delegate.h"
#include "xml_parameter_reader.h"
#include "xml_parameter_writer.h"


namespace dealii
{
  namespace ParameterGui
  {
    MainWindow::MainWindow(const QString  &filename)
    {
      QString  settings_file = QDir::currentPath() + "/settings.ini";		// a file for user settings

      gui_settings = new QSettings (settings_file, QSettings::IniFormat);	// load settings
										// Up to now, we do not read any settings,
										// but this can be used in the future for customizing the GUI.

      tree_widget = new QTreeWidget;						// tree for showing XML tags

										// Setup the tree and the window first:
      tree_widget->header()->setResizeMode(QHeaderView::ResizeToContents);	// behavior of the header sections:
										// "Interactive: User can resize sections"
										// "Fixed: User cannot resize sections"
										// "Stretch: Qt will automatically resize sections to fill available space"
										// "ResizeToContents: Qt will automatically resize sections to optimal size"
      tree_widget->setHeaderLabels(QStringList() << tr("(Sub)Sections/Parameters")
                                                 << tr("Value"));
      tree_widget->setMouseTracking(true);					// enables mouse events e.g. showing ToolTips
										// and documentation in the StatusLine
      tree_widget->setEditTriggers(QAbstractItemView::DoubleClicked|
                                   QAbstractItemView::SelectedClicked|
                                   QAbstractItemView::EditKeyPressed);
										// set which actions will initiate item editing: Editing starts when:
										// DoubleClicked: an item is double clicked
										// SelectedClicked: clicking on an already selected item
										// EditKeyPressed: the platform edit key has been pressed over an item
										// AnyKeyPressed: any key is pressed over an item

      tree_widget->setItemDelegate(new ParameterDelegate(1));			// set the delegate for editing items
      setCentralWidget(tree_widget);
										// connect: if the tree changes, the window will know
      connect(tree_widget, SIGNAL(itemChanged(QTreeWidgetItem*, int)), this, SLOT(tree_was_modified()));

      create_actions();								// create window actions as "Open",...
      create_menus();								// and menus
      statusBar()->showMessage(tr("Ready, start editing by double-clicking or hitting F2!"));
      setWindowTitle(tr("[*]parameterGUI"));					// set window title

      resize(800, 600);								// set window height and width

      if (filename.size() > 3)							// if there is a file_name, try to load the file.
        load_file(filename);							// a vliad file has the xml extension, so we require size() > 3
    }



    void MainWindow::open()
    {
      if (maybe_save())								// check, if the content was modified
        {
          QString  file_name =							// open a file dialog
                     QFileDialog::getOpenFileName(this, tr("Open XML Parameter File"),
                                                  QDir::currentPath(),
                                                  tr("XML Files (*.xml)"));
          if (!file_name.isEmpty())						// if a file was selected,
            load_file(file_name);						// load the content
        };
    }



    bool MainWindow::save()
    {
      if (current_file.isEmpty())						// if there is no file
        return save_as();							// to save changes, open a dialog
      else
        return save_file(current_file);						// otherwise save
    }



    bool MainWindow::save_as()
    {
      QString  file_name =							// open a file dialog
                 QFileDialog::getSaveFileName(this, tr("Save XML Parameter File"),
                                              QDir::currentPath(),
                                              tr("XML Files (*.xml)"));

      if (file_name.isEmpty())							// if no file was selected
        return false;								// return false
      else
        return save_file(file_name);						// otherwise save content to file
    }



    void MainWindow::about()
    {
#ifdef Q_WS_MAC
      static QPointer<QMessageBox> old_msg_box;

      if (old_msg_box)
        {
          old_msg_box->show();
          old_msg_box->raise();
          old_msg_box->activateWindow();
          return;
        };
#endif

      QString title = "About parameterGUI";

      QString trAboutparameterGUIcaption;
      trAboutparameterGUIcaption = QMessageBox::tr(
        "<h3>parameterGUI: A GraphicalUserInterface for parameter handling in deal.II</h3>"
        "<p>This program uses Qt version %1.</p>"
        ).arg(QLatin1String(QT_VERSION_STR));

      QString trAboutparameterGUItext;
      trAboutparameterGUItext = QMessageBox::tr(
        "<p>The parameterGUI is a graphical user interface for editing XML parameter files "
        "created by the ParameterHandler class of deal.II. Please see "
        "<a href=\"http://www.dealii.org/7.0.0/doxygen/deal.II/classParameterHandler.html\">dealii.org/doc</a> for more information. "
        "The parameterGUI parses XML files into a tree structure and provides "
        " special editors for different types of parameters.</p>"

        "<p><b>Editing parameter values:</b><br>"
        "Parameters can be edited by (double-)clicking on the value or "
        "by pressing the platform edit key (F2 on Linux) over an parameter item.</p>"

        "<p><b>Editors for parameter values:</b>"
        " <ul>"
        "  <li>Integer- and Double-type parameters: SpinBox</li>"
        "  <li>Booleans: ComboBox</li>"
        "  <li>Selection: ComboBox</li>"
        "  <li>File- and DirectoryName parameters: BrowseLineEditor</li>"
        "  <li>Anything|MultipleSelection|List: LineEditor</li>"
        " </ul>"
        "</p>"

        "<p>Please see <a href=\"http://www.dealii.org\">dealii.org</a> for more information</p>"
        "<p><b>Authors:</b><br> "
        "Martin Steigemann,  <a href=\"mailto:martin.steigemann@mathematik.uni-kassel.de\">martin.steigemann@mathematik.uni-kassel.de</a><br>"
        "Wolfgang Bangerth,  <a href=\"mailto:bangerth@math.tamu.edu\">bangerth@math.tamu.edu</a></p>"
        );

      QMessageBox *msg_box = new QMessageBox;
      msg_box->setAttribute(Qt::WA_DeleteOnClose);
      msg_box->setWindowTitle(title);
      msg_box->setText(trAboutparameterGUIcaption);
      msg_box->setInformativeText(trAboutparameterGUItext);

      QPixmap pm(QLatin1String(":/images/logo_dealii_gui_128.png"));

      if (!pm.isNull())
        msg_box->setIconPixmap(pm);

#ifdef Q_WS_MAC
      old_msg_box = msg_box;
      msg_box->show();
#else
      msg_box->exec();
#endif
    }



    void MainWindow::tree_was_modified()
    {
      setWindowModified(true);							// store, that the window was modified
										// this is a function from the QMainWindow class
										// and we use the windowModified mechanism to show a "*"
										// in the window title, if content was modified
    }



    void MainWindow::show_message ()
    {
      QString title = "parameterGUI";

      info_message = new InfoMessage(this);

      info_message->setWindowTitle(title);
      info_message->setInfoMessage(tr("Start Editing by double-clicking on the parameter value or"
                                      " by hitting the platform edit key. For example, on Linux this is the F2-key!"));
      info_message->showMessage();
    }



    void MainWindow::closeEvent(QCloseEvent *event)
    {
      if (maybe_save())								// reimplement the closeEvent from the QMainWindow class
        event->accept();							// check, if we have to save modified content,
      else									// if content was saved, accept the event,
        event->ignore();							// otherwise ignore it
    }



    void MainWindow::create_actions()
    {
      QStyle * style = tree_widget->style();

      open_act = new QAction(tr("&Open..."), this);				// create actions
      open_act->setIcon(style->standardPixmap(QStyle::SP_DialogOpenButton));    // and set icons
      open_act->setShortcut(Qt::CTRL + Qt::Key_O);				// set a short cut
      open_act->setStatusTip(tr("Open a XML file"));				// set a status tip
      connect(open_act, SIGNAL(triggered()), this, SLOT(open()));		// and connect

      save_act = new QAction(tr("&Save ..."), this);
      save_act->setIcon(style->standardPixmap(QStyle::SP_DialogSaveButton));
      save_act->setShortcut(Qt::CTRL + Qt::Key_S);
      save_act->setStatusTip(tr("Save the current XML file"));
      connect(save_act, SIGNAL(triggered()), this, SLOT(save()));

      save_as_act = new QAction(tr("&Save As..."), this);
      save_as_act->setIcon(style->standardPixmap(QStyle::SP_DialogSaveButton));
      save_as_act->setShortcut(Qt::CTRL + Qt::SHIFT + Qt::Key_Q);
      save_as_act->setStatusTip(tr("Save the current XML file as"));
      connect(save_as_act, SIGNAL(triggered()), this, SLOT(save_as()));

      exit_act = new QAction(tr("E&xit"), this);
      exit_act->setIcon(style->standardPixmap(QStyle::SP_DialogCloseButton));
      exit_act->setShortcut(Qt::CTRL + Qt::Key_Q);
      exit_act->setStatusTip(tr("Exit the parameterGUI application"));
      connect(exit_act, SIGNAL(triggered()), this, SLOT(close()));

      about_act = new QAction(tr("&About"), this);
      about_act->setIcon(style->standardPixmap(QStyle::SP_FileDialogInfoView));
      about_act->setStatusTip(tr("Show the parameterGUI About box"));
      connect(about_act, SIGNAL(triggered()), this, SLOT(about()));

      about_qt_act = new QAction(tr("About &Qt"), this);
      about_qt_act->setStatusTip(tr("Show the Qt library's About box"));
      connect(about_qt_act, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
    }



    void MainWindow::create_menus()
    {
        file_menu = menuBar()->addMenu(tr("&File"));				// create a file menu
        file_menu->addAction(open_act);						// and add actions
        file_menu->addAction(save_act);
        file_menu->addAction(save_as_act);
        file_menu->addAction(exit_act);

        menuBar()->addSeparator();

        help_menu = menuBar()->addMenu(tr("&Help"));				// create a help menu
        help_menu->addAction(about_act);
        help_menu->addAction(about_qt_act);
    }



    bool MainWindow::maybe_save()
    {
      if (isWindowModified())							// if content was modified
        {
          QMessageBox::StandardButton ret;					// ask, if content should be saved
          ret = QMessageBox::warning(this, tr("parameterGUI"),
                                     tr("The content has been modified.\n"
                                        "Do you want to save your changes?"),
                  QMessageBox::Save | QMessageBox::Discard |QMessageBox::Cancel);

          if (ret == QMessageBox::Save)
            return save();
          else if (ret == QMessageBox::Cancel)
            return false;
        };

      return true;
    }



    bool MainWindow::save_file(const QString &filename)
    {
      QFile  file(filename);

      if (!file.open(QFile::WriteOnly | QFile::Text))				// open a file dialog
        {
          QMessageBox::warning(this, tr("parameterGUI"),
                                     tr("Cannot write file %1:\n%2.")
                                     .arg(filename)
                                     .arg(file.errorString()));
          return false;
        };

      XMLParameterWriter  xml_writer(tree_widget);				// create a xml_writer

      if (!xml_writer.write_xml_file(&file))					// and read the xml file
        return false;

      statusBar()->showMessage(tr("File saved"), 2000);				// if we succeed, show a message
      set_current_file(filename);						// and reset the window

      return true;
    }



    void MainWindow::load_file(const QString &filename)
    {
      QFile  file(filename);

      if (!file.open(QFile::ReadOnly | QFile::Text))				// open the file
        {
          QMessageBox::warning(this, tr("parameterGUI"),
                                     tr("Cannot read file %1:\n%2.")
                                     .arg(filename)
                                     .arg(file.errorString()));
          return;
        };

      tree_widget->clear();							// clear the tree

      XMLParameterReader  xml_reader(tree_widget);				// and read the xml file

      if (!xml_reader.read_xml_file(&file))
        {
          QMessageBox::warning(this, tr("parameterGUI"),
                                     tr("Parse error in file %1:\n\n%2")
                                     .arg(filename)
                                     .arg(xml_reader.error_string()));
        }
      else
        {
          statusBar()->showMessage(tr("File loaded - Start editing by double-clicking or hitting F2"), 25000);
          set_current_file(filename);						// show a message and set current file

          show_message ();							// show some informations how values can be edited
        };
    }



    void MainWindow::set_current_file(const QString  &filename)
    {
										// We use the windowModified mechanism from the
										// QMainWindow class to indicate in the window title,
										// if the content was modified.
										// If there is "[*]" in the window title, a * will
										// added automatically at this position, if the
										// window was modified.
										// We set the window title to
										// file_name[*] - XMLParameterHandler

      current_file = filename;							// set the (global) current file to file_name

      std::string win_title = (filename.toStdString());				// and create the window title,

      if (current_file.isEmpty())						// if file_name is empty
        win_title = "[*]parameterGUI";						// set the title to our application name,
      else
        win_title += "[*] - parameterGUI";					// if there is a file_name, add the
										// the file_name and a minus to the title

      setWindowTitle(tr(win_title.c_str()));					// set the window title
      setWindowModified(false);							// and reset window modified
    }
  }
}
