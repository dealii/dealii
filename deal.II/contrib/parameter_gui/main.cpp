
#include <QApplication>
#include <QSplashScreen>
#include <QTimer>

#include "mainwindow.h"


/**
 * Main function for the parameterGUI.
 * The parameterGUI is a graphical user interface for editing parameter files based on the XML format,
 * created by the ParameterHandler::print_parameters() function with ParameterHandler::XML as second argument.
 *
 *   @image html logo_dealii_gui.png
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
 */
int main(int argc, char *argv[])
{
  Q_INIT_RESOURCE(application);					// init resources such as icons or graphics

  QApplication app(argc, argv);

  QSplashScreen * splash = new QSplashScreen;			// setup a splash screen
  splash->setPixmap(QPixmap(":/images/logo_dealii_gui.png"));
  splash->show();

  QTimer::singleShot(3000, splash, SLOT(close()));		// and close it after 3000 ms

  app.setApplicationName("parameterGUI for deal.II");		// setup the application name

  MainWindow * main_win = new MainWindow (argv[1]);		// give command line arguments to main_win
								// if a parameter file is specified at the
								// command line, give it to the MainWindow.

  QTimer::singleShot(1500, main_win, SLOT(show()));		// show the main window with a short delay
								// so we can see the splash screen
  return app.exec();
}
