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


#ifndef BROWSELINEEDIT_H
#define BROWSELINEEDIT_H

#include <QFrame>
#include <QLineEdit>
#include <QFileDialog>
#include <QPushButton>


namespace dealii
{
/*! @addtogroup ParameterGui
 *@{
 */
  namespace ParameterGui
  {
/**
 * The BrowseLineEdit class provides a special line editor for the parameterGUI.
 * While editing file- or directory names it is much more easier to have a file-dialog
 * and just click on existing files or directories. This editor provides a simple QLineEditor
 * and a browse-button which opens a file- or a directory dialog. Clicking on existing files or directories
 * copies the path to the line editor. Depending on the <tt>BrowseType</tt> given in the constructor
 * the browse button opens a <tt>file</tt> or a <tt>directory</tt> dialog.
 *
 * @note This class is used in the graphical user interface for the @ref ParameterHandler class.
 *       It is not compiled into the deal.II libraries and can not be used by applications using deal.II.
 *
 * @ingroup ParameterGui
 * @author Martin Steigemann, Wolfgang Bangerth, 2010
 */
    class BrowseLineEdit : public QFrame
    {
      Q_OBJECT

      public:
				     /**
				      * The browse button opens a <tt>file</tt> or
				      * a <tt>directory</tt> dialog. This can be specified
				      * in the constructor by setting this flag <tt>BrowseType</tt>.
				      */
        enum BrowseType {file = 0, directory = 1};
				     /**
				      * Constructor. The type of the browse dialog can be specified
				      * by the flag <tt>BrowseType</tt>, the default is <tt>file</tt>.
				      */
        BrowseLineEdit (const BrowseType  type = file,
                        QWidget          *parent = 0);

				     /**
				      * Reimplemented from the QWidget class.
				      * Returns the size of the editor.
				      */
        QSize  sizeHint() const;
				     /**
				      * Reimplemented from the QWidget class.
				      */
        QSize  minimumSizeHint() const;
				     /**
				      * Returns the text of the line editor.
				      */
        QString  text() const;
				     /**
				      * This pattern stores the type of the browse dialog.
				      */
        BrowseType  browse_type;

      public slots:
				     /**
				      * A <tt>slot</tt> to set @p str as text of the line editor.
				      */
        void setText(const QString &str);

      signals:
				     /**
				      * This <tt>signal</tt> will be emitted, if editing is finished.
				      */
        void editingFinished();

      private slots:
				     /**
				      * This <tt>slot</tt> should be always called, if editing is finished.
				      */
        void editing_finished();
				     /**
				      * This function opens a file- or a directory dialog as specified in the
				      * constructor.
				      */
        void browse();

      private:
				     /**
				      * The line editor.
				      */
        QLineEdit * line_editor;
				     /**
				      * The browse button.
				      */
        QPushButton * browse_button;
    };
  }
/**@}*/
}


#endif
