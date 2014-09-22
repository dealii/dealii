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

#include "info_message.h"


namespace dealii
{
  namespace ParameterGui
  {
    InfoMessage::InfoMessage(QWidget *parent)
               : QDialog(parent, 0)
    {
      show_again = true;						// this variable stores, if the
									// the info message should be shown again
      QGridLayout * grid = new QGridLayout(this);

      icon = new QLabel(this);						// set an icon
#ifndef QT_NO_MESSAGEBOX
      icon->setPixmap(QMessageBox::standardIcon(QMessageBox::Information));
      icon->setAlignment(Qt::AlignHCenter | Qt::AlignTop);
#endif
      grid->addWidget(icon, 0, 0, Qt::AlignTop);			// add the icon in the upper left corner

      message = new QTextEdit(this);					// set the new message
      message->setReadOnly(true);
      grid->addWidget(message, 0, 1);					// and add the message on the right

      again = new QCheckBox(this);					// add a check box
      again->setChecked(true);
      again->setText(QErrorMessage::tr("&Show this message again"));
      grid->addWidget(again, 1, 1, Qt::AlignTop);

      ok = new QPushButton(this);					// and finally a OK button
      ok->setText(QErrorMessage::tr("&OK"));
#ifdef QT_SOFTKEYS_ENABLED
      ok_action = new QAction(ok);					// define the action for the button
      ok_action->setSoftKeyRole(QAction::PositiveSoftKey);
      ok_action->setText(ok->text());
      connect(ok_action, SIGNAL(triggered()), this, SLOT(accept()));
      addAction(ok_action);
#endif
      connect(ok, SIGNAL(clicked()), this, SLOT(accept()));
      ok->setFocus();							// aand set the focus on the button
      grid->addWidget(ok, 2, 0, 1, 2, Qt::AlignCenter);

      grid->setColumnStretch(1, 42);
      grid->setRowStretch(0, 42);
									// load settings from an ini-file
      QString  settings_file = QDir::currentPath() + "/settings.ini";

      settings = new QSettings (settings_file, QSettings::IniFormat);

      settings->beginGroup("infoMessage");				// we store settings of this class in the
      show_again = settings->value("showInformation", true).toBool();	//group infoMessage
      settings->endGroup();
    }



    void InfoMessage::setInfoMessage(const QString &message)
    {
      this->message->setText(message);					// set the message
    }



    void InfoMessage::showMessage()
    {
      if (show_again)							// and show the message
        show();
    }



    void InfoMessage::done(int r)
    {
      if(!again->isChecked())						// if the box is not checked,
        {								// store this to settings
          settings->beginGroup("infoMessage");
          settings->setValue("showInformation", false);
          settings->endGroup();
        };

      QDialog::done(r);
    }
  }
}

