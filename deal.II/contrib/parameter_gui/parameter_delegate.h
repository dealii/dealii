
#ifndef PARAMETERDELEGATE_H
#define PARAMETERDELEGATE_H

#include <QItemDelegate>
#include <QModelIndex>
#include <QObject>
#include <QLineEdit>
#include <QComboBox>
#include <QFileDialog>

#include <limits>
#include <algorithm>

#include "browse_lineedit.h"


/**
 * The ParameterDelegate class implements special delegates for the QTreeWidget class used in the parameterGUI.
 * The QTreeWidget class provides some different standard delegates for editing parameters shown in the
 * tree structure. The ParameterDelegate class provides special editors for the different types of parameters defined in
 * the ParameterHandler class. For all parameter types based on strings as "Anything", "MultipleSelection" and
 * "List" a simple line editor will be shown up. In the case of integer and double type parameters the editor is a spin box and for
 * "Selection" type parameters a combo box will be shown up. For parameters of type "FileName" and "DirectoryName"
 * the delegate shows a @ref BrowseLineEdit editor. The column of the tree structure with the parameter values has to be set
 * in the constructor.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
 */
class ParameterDelegate : public QItemDelegate
{
  Q_OBJECT

  public:
				     /**
				      * Constructor, the <tt>value_column</tt> specifies the column
				      * of the parameter tree this delegate will be used on.
				      */
    ParameterDelegate (const int value_column, QObject *parent = 0);
				     /**
				      * This function creates the appropriate editor for the parameter
				      * based on the <tt>index</tt>.
				      */
    QWidget * createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                           const QModelIndex &index) const;
				     /**
				      * Reimplemented from QItemDelegate.
				      */
    QSize sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const;
				     /**
				      * Reimplemented from QItemDelegate.
				      */
    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;
				     /**
				      * Reimplemented from QItemDelegate.
				      */
    void setEditorData(QWidget *editor, const QModelIndex &index) const;
				     /**
				      * Reimplemented from QItemDelegate.
				      */
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const;

  private slots:
				     /**
				      * Reimplemented from QItemDelegate.
				      */
    void commit_and_close_editor();

  private:
				     /**
				      * The column this delegate will be used on.
				      */
    int value_column;
				     /**
				      * For parameters of type <tt>double</tt> a spin box
				      * will be shown as editor. Any click on the spin box
				      * will change the value about <tt>double_steps</tt>.
				      */
    double  double_steps;
				     /**
				      * For parameters of type <tt>integer</tt> a spin box
				      * will be shown as editor. Any click on the spin box
				      * will change the value about <tt>int_steps</tt>.
				      */
    unsigned int  int_steps;
				     /**
				      * For parameters of type <tt>double</tt> a spin box
				      * will be shown as editor. The spin box will show
				      * parameters with a precision of <tt>double_decimals</tt>.
				      */
    unsigned int  double_decimals;
};

#endif
