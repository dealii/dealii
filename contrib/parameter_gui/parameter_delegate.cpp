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

#include "parameter_delegate.h"

#include <limits>


namespace dealii
{
  namespace ParameterGui
  {
    ParameterDelegate::ParameterDelegate(const int value_column, QObject *parent)
                     : QItemDelegate(parent)
    {
      this->value_column = value_column;

      double_steps    = 0.1;			// any click in the editor will increase or decrease the value about double_steps
      double_decimals = 14;			// number of decimals shown in the editor

      int_steps = 1;				// step value for increasing or decrasing integers
    }



    QSize ParameterDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const
    {
      if (index.column() == value_column)
        {
          return QSize(400,30);		// we increase the height of all lines to show editors

/*
      QString pattern_description = index.data(Qt::StatusTipRole).toString();	// load pattern description
										// stored in the StatusLine
      QRegExp  rx_string("\\b(FileName|DirectoryName)\\b");

      if (rx_string.indexIn (pattern_description) != -1)
        {
          return QSize(400,35);					// we increase the height for FileName and
        }							// DirectoryName to show a "browse" button
      else
        return QItemDelegate::sizeHint(option, index);
*/

        }
      else
        return QItemDelegate::sizeHint(option, index);
    }



    void ParameterDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
    {
      if (index.column() == value_column)
        {
          QString pattern_description = index.data(Qt::StatusTipRole).toString();	// load pattern description
											// stored in the StatusLine
          QRegExp  rx_string("\\b(FileName|DirectoryName)\\b");				// if the type is Filename
											// or DirectoryName
          if (rx_string.indexIn (pattern_description) != -1)
            {
              QString value = index.model()->data(index, Qt::DisplayRole).toString();	// take the value

              QStyleOptionViewItem my_option = option;					// load options
              my_option.displayAlignment = Qt::AlignLeft | Qt::AlignVCenter;

              drawDisplay(painter, my_option, my_option.rect, value);			// print the text in the display
              drawFocus(painter, my_option, my_option.rect);				// if the line has the
											// focus, print a rectangle
            }
          else
            QItemDelegate::paint(painter, option, index);				// for all other types use
											// the standard delegate
        }
      else
        QItemDelegate::paint(painter, option, index);
    }



    QWidget *ParameterDelegate::createEditor(QWidget *parent,
                                             const QStyleOptionViewItem &option,
                                             const QModelIndex &index) const
    {
      if (index.column() == value_column)
        {
          QString pattern_description = index.data(Qt::StatusTipRole).toString();	// load pattern description
											// stored in the StatusLine
          QRegExp  rx_string("\\b(Anything|MultipleSelection|List|Map)\\b"),
                   rx_filename("\\b(FileName)\\b"),	
                   rx_dirname("\\b(DirectoryName)\\b"),	
                   rx_integer("\\b(Integer)\\b"),
                   rx_double("\\b(Double|Float|Floating)\\b"),
                   rx_selection("\\b(Selection)\\b"),
                   rx_bool("\\b(Bool)\\b");

          if (rx_string.indexIn (pattern_description) != -1)				// if the type is "Anything"
            {
              QLineEdit * line_editor = new QLineEdit(parent);				// choose a LineEditor

              connect(line_editor, SIGNAL(editingFinished()),				// and connect editors signal
                      this, SLOT(commit_and_close_editor()));				// to the closer function

              return line_editor;
            }
          else if (rx_filename.indexIn (pattern_description) != -1)			// if the type is "FileName"
            {
              BrowseLineEdit * filename_editor =					// choose a BrowseLineEditor
                                 new BrowseLineEdit(BrowseLineEdit::file, parent);

              connect(filename_editor, SIGNAL(editingFinished()),
                      this, SLOT(commit_and_close_editor()));

              return filename_editor;
            }
          else if (rx_dirname.indexIn (pattern_description) != -1)			// if the type is "DirectoryName"
            {
              BrowseLineEdit * dirname_editor =						// choose a BrowseLineEditor
                                 new BrowseLineEdit(BrowseLineEdit::directory, parent);

              connect(dirname_editor, SIGNAL(editingFinished()),
                      this, SLOT(commit_and_close_editor()));

              return dirname_editor;
            }
          else if (rx_integer.indexIn (pattern_description) != -1)		// if the tpye is "Integer"
            {
              QSpinBox * spin_box = new QSpinBox(parent);			// choose a spin box

              const int min_int_value = std::numeric_limits<int>::min();
              const int max_int_value = std::numeric_limits<int>::max();

              spin_box->setMaximum(max_int_value);				// set max and min from the limits.h class
              spin_box->setMinimum(min_int_value);
              spin_box->setSingleStep(int_steps);				// and every klick is a SingleStep

              connect(spin_box, SIGNAL(editingFinished()),			// connect editors signal to the closer function
                      this, SLOT(commit_and_close_editor()));

              return spin_box;
            }
          else if (rx_double.indexIn (pattern_description) != -1)		// the same with "Double"
            {
              QDoubleSpinBox * double_spin_box = new QDoubleSpinBox(parent);	// choose a spin box

              const double min_double_value = -std::numeric_limits<double>::max();
              const double max_double_value = std::numeric_limits<double>::max();

              double_spin_box->setMaximum(max_double_value);		// set max and min from the limits.h class
              double_spin_box->setMinimum(min_double_value);
              double_spin_box->setDecimals(double_decimals);		// show "double_decimals" decimals
              double_spin_box->setSingleStep(double_steps);		// and every klick is a SingleStep

              connect(double_spin_box, SIGNAL(editingFinished()),		// connect editors signal to the closer function
                      this, SLOT(commit_and_close_editor()));

              return double_spin_box;
            }
          else if (rx_selection.indexIn (pattern_description) != -1)		// and selections
            {
              QComboBox * combo_box = new QComboBox(parent);			// we assume, that pattern_desctiption is of the form
										// "Type: [Selection item1|item2| ....|item ]    "
              std::vector<std::string> choices;					// list with the different items
              std::string  tmp(pattern_description.toStdString());

              if (tmp.find("[") != std::string::npos)				// delete all char before [
                tmp.erase (0, tmp.find("[")+1);

              if (tmp.find("]") != std::string::npos)				// delete all char after ]
                tmp.erase (tmp.find("]"),tmp.length());

              if (tmp.find(" ") != std::string::npos)				// delete all char before " "
                tmp.erase (0, tmp.find(" ")+1);

              while (tmp.find('|') != std::string::npos)			// extract items
                {
                  choices.push_back(std::string(tmp, 0, tmp.find('|')));
                  tmp.erase (0, tmp.find('|')+1);
                };

              if (tmp.find(" ") != std::string::npos)				// delete " "
                tmp.erase (tmp.find(" "));

              choices.push_back(tmp);						// add last item

              for (unsigned int i=0; i<choices.size(); ++i)			// add items to the combo box
                combo_box->addItem (tr(choices[i].c_str()), tr(choices[i].c_str()));

              combo_box->setEditable(false);

              connect(combo_box, SIGNAL(currentIndexChanged(int)),		// connect editors signal to the closer function
                      this, SLOT(commit_and_close_editor()));

              return combo_box;
           }
          else if (rx_bool.indexIn (pattern_description) != -1)			// and booleans
            {
              QComboBox * combo_box = new QComboBox(parent);

              std::vector<std::string> choices;					// list with the different items
              choices.push_back(std::string("true"));				// add true
              choices.push_back(std::string("false"));				// and false

              for (unsigned int i=0; i<choices.size(); ++i)			// add items to the combo box
                combo_box->addItem (tr(choices[i].c_str()), tr(choices[i].c_str()));

              combo_box->setEditable(false);

              connect(combo_box, SIGNAL(currentIndexChanged(int)),		// connect editors signal to the closer function
                      this, SLOT(commit_and_close_editor()));

              return combo_box;
            }
          else
            {
              return QItemDelegate::createEditor(parent, option, index);
            };
        };

      return 0;				// if it is not the column "parameter values", do nothing
    }



    void ParameterDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
    {
      if (index.column() == value_column)
        {
          QString pattern_description = index.data(Qt::StatusTipRole).toString();	// load pattern description
											// stored in the StatusLine
          QRegExp  rx_filename("\\b(FileName)\\b"),
                   rx_dirname("\\b(DirectoryName)\\b"),
                   rx_selection("\\b(Selection)\\b");

          if (rx_filename.indexIn (pattern_description) != -1)				// if the type is "FileName"
            {
              QString  file_name = index.data(Qt::DisplayRole).toString();

              BrowseLineEdit *filename_editor = qobject_cast<BrowseLineEdit *>(editor);	// set the text of the editor
              filename_editor->setText(file_name);
            }
          else if (rx_dirname.indexIn (pattern_description) != -1)			// if the type is "DirectoryName"
            {
              QString  dir_name = index.data(Qt::DisplayRole).toString();

              BrowseLineEdit *dirname_editor = qobject_cast<BrowseLineEdit *>(editor);	// set the text of the editor
              dirname_editor->setText(dir_name);
            }
          else if (rx_selection.indexIn (pattern_description) != -1)			// if we have a combo box,
            {
              QRegExp  rx(index.data(Qt::DisplayRole).toString());

              QComboBox * combo_box = qobject_cast<QComboBox *>(editor);

              for (int i=0; i<combo_box->count(); ++i)					// we look, which index
                if (rx.exactMatch(combo_box->itemText(i)))				// the data has and set
                  combo_box->setCurrentIndex(i);					// it to the combo_box
            }
          else
            QItemDelegate::setEditorData(editor, index);				// if it is not FileName,
											// DirectoryName or Selection
											// use the standard delegate
        };
    }



    void ParameterDelegate::commit_and_close_editor()
    {
      QWidget * editor = qobject_cast<QWidget *>(sender());
      emit commitData(editor);
      emit closeEditor(editor);
    }



    void ParameterDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                         const QModelIndex &index) const
    {
      if (index.column() == value_column)
        {
          QString pattern_description = index.data(Qt::StatusTipRole).toString();	// load pattern description
											// stored in the StatusLine

          QRegExp  rx_filename("\\b(FileName)\\b"),
                   rx_dirname("\\b(DirectoryName)\\b"),
                   rx_selection("\\b(Selection)\\b");

          if (rx_filename.indexIn (pattern_description) != -1)				// if the type is "FileName"
            {
              BrowseLineEdit * filename_editor = qobject_cast<BrowseLineEdit *>(editor);	// set the text from the editor
              QString value = filename_editor->text();
              model->setData(index, value);
            }
          else if (rx_dirname.indexIn (pattern_description) != -1)			// if the type is "DirectoryName"
            {
              BrowseLineEdit * dirname_editor = qobject_cast<BrowseLineEdit *>(editor);	// set the text from the editor
              QString value = dirname_editor->text();
              model->setData(index, value);
            }
          else if (rx_selection.indexIn (pattern_description) != -1)			// if the type is "Selection"
            {
              QComboBox * combo_box = qobject_cast<QComboBox *>(editor);		// set the text from the editor
              QString value = combo_box->currentText();
              model->setData(index, value);
            }
          else
            QItemDelegate::setModelData(editor, model, index);				// if it is not FileName or DirectoryName,
											// use the standard delegate
        };
    }
  }
}

