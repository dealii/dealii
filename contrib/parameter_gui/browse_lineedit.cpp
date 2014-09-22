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

#include "browse_lineedit.h"


namespace dealii
{
  namespace ParameterGui
  {
    BrowseLineEdit::BrowseLineEdit(const BrowseType  type, QWidget *parent)
                  : QFrame(parent, 0),
                    browse_type(type)
    {
      line_editor = new QLineEdit;
      connect(line_editor, SIGNAL(editingFinished()), this, SLOT(editing_finished()));

      browse_button = new QPushButton("&Browse...");
      connect(browse_button, SIGNAL(clicked()), this, SLOT(browse()));

      setFocusPolicy (Qt::StrongFocus);

      QHBoxLayout *layout = new QHBoxLayout;

      layout->addWidget(line_editor);
      layout->addWidget(browse_button);
      setLayout(layout);

      setAutoFillBackground(true);
      setBackgroundRole(QPalette::Highlight);
    }




    QSize BrowseLineEdit::sizeHint() const
    {
      QSize  size_line_editor   = line_editor->sizeHint(),
             size_browse_button = browse_button->sizeHint();

      int w = size_line_editor.rwidth() + size_browse_button.rwidth(),
          h = qMax(size_line_editor.rheight(), size_browse_button.rheight());

      return QSize (w, h);
    }



    QSize BrowseLineEdit::minimumSizeHint() const
    {
      QSize  size_line_editor   = line_editor->minimumSizeHint(),
             size_browse_button = browse_button->minimumSizeHint();

      int w = size_line_editor.rwidth() + size_browse_button.rwidth(),
          h = qMax(size_line_editor.rheight(), size_browse_button.rheight());

      return QSize (w, h);
    }



    QString BrowseLineEdit::text() const
    {
      return line_editor->text();
    }



    void BrowseLineEdit::setText(const QString &str)
    {
      line_editor->setText(str);
    }



    void BrowseLineEdit::editing_finished()
    {
      emit editingFinished();
    }



    void BrowseLineEdit::browse()
    {
      QString  name = "";

      switch (browse_type)
        {
          case file:
            {
              name = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                  QDir::currentPath(),
                                                  tr("All Files (*.*)"));
              break;
            };

          case directory:
            {
              name = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                       QDir::homePath(),
                                                       QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
              break;
            };
        };

      if (!name.isEmpty() && !name.isNull())
        line_editor->setText(name);
    }
  }
}

