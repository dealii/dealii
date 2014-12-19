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


#ifndef INFOMESSAGE_H
#define INFOMESSAGE_H

#include <QDialog>
#include <QSettings>
#include <QCheckBox>
#include <QTextEdit>
#include <QLabel>


namespace dealii
{
/*! @addtogroup ParameterGui
 *@{
 */
  namespace ParameterGui
  {
/**
 * The InfoMessage class implements a special info message for the parameterGUI.
 * Besides showing a info message itself, the dialog shows a checkbox &quot;Show this message again&quot;.
 * If the user unchecks this box, this is stored in the &quot;settings.ini&quot; file and will be reloaded
 * the next time the user opens the parameterGUI. The intention of such a info message is the following.
 * The user should have some information on how using the GUI &quot;at hand&quot;
 * such as &quot;how to edit parameter values&quot; for example. But after reading this message, the user knows
 * it and the message should not appear permanently.
 *
 * @note This class is used in the graphical user interface for the @ref ParameterHandler class.
 *       It is not compiled into the deal.II libraries and can not be used by applications using deal.II.
 *
 * @ingroup ParameterGui
 * @author Martin Steigemann, Wolfgang Bangerth, 2010
 */
    class InfoMessage : public QDialog
    {
      Q_OBJECT

      public:
				     /**
				      * Constructor
				      */
        InfoMessage (QWidget *parent = 0);
				     /**
				      * With this function the @p message which will be shown in the
				      * dialog can be set.
				      */
        void setInfoMessage(const QString &message);

      public slots:
				     /**
				      * Show the dialog with the <tt>message</tt>.
				      */
        void showMessage();

      protected:
				     /**
				      * Reimplemented from QDialog.
				      */
        void done(int r);

      private:
				     /**
				      * This variable stores, if the <tt>message</tt> should be shown again the next time.
				      */
        bool show_again;
				     /**
				      * The <tt>Ok</tt> button.
				      */
        QPushButton * ok;
				     /**
				      * The checkbox<tt>Show this message again</tt>.
				      */
        QCheckBox * again;
				     /**
				      * The <tt>message</tt> editor.
				      */
        QTextEdit * message;
				     /**
				      * An <tt>icon</tt> for the dialog.
				      */
        QLabel * icon;
#ifdef QT_SOFTKEYS_ENABLED
				     /**
				      * A action for pressing the <tt>Ok</tt> button.
				      */
        QAction * ok_action;
#endif
				     /**
				      * An object for storing <tt>settings</tt> in a file.
				      */
        QSettings * settings;
    };
  }
/**@}*/
}


#endif
