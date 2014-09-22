// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by Martin Steigemann and Wolfgang Bangerth
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


#include <QApplication>
#include <QSplashScreen>
#include <QTimer>

#include "mainwindow.h"

/*! @addtogroup ParameterGui
 *@{
 */

/**
 * Main function for the parameterGUI.
 * The parameterGUI is a graphical user interface for editing parameter files based on the XML format,
 * created by the ParameterHandler::print_parameters() function with ParameterHandler::XML as second argument.
 *
 * @image html logo_dealii_gui.png
 *
 * @note This class is used in the graphical user interface for the @ref ParameterHandler class.
 *       It is not compiled into the deal.II libraries and can not be used by applications using deal.II.
 *
 *
 * <p>This program uses Qt version > 4.3. Qt is licensed under the GNU General Public License
 * version 3.0. Please see http://qt.nokia.com/products/licensing for an overview of Qt licensing.
 * Copyright (C) 2010 Nokia Corporation and/or its subsidiary(-ies). Qt is a Nokia product.
 * See http://qt.nokia.com/ for more information.</p>
 *
 *
 * @ingroup ParameterGui
 * @author Martin Steigemann, Wolfgang Bangerth, 2010
 */
int main(int argc, char *argv[])
{
  Q_INIT_RESOURCE(application);						// init resources such as icons or graphics

  QApplication app(argc, argv);

  QSplashScreen * splash = new QSplashScreen;				// setup a splash screen
  splash->setPixmap(QPixmap(":/images/logo_dealii_gui.png"));
  splash->show();

  QTimer::singleShot(3000, splash, SLOT(close()));			// and close it after 3000 ms

  app.setApplicationName("parameterGUI for deal.II");			// setup the application name

  dealii::ParameterGui::MainWindow * main_win =
    new dealii::ParameterGui::MainWindow (argv[1]);			// give command line arguments to main_win
									// if a parameter file is specified at the
									// command line, give it to the MainWindow.

  QTimer::singleShot(1500, main_win, SLOT(show()));			// show the main window with a short delay
									// so we can see the splash screen
  return app.exec();
}
/**@}*/

