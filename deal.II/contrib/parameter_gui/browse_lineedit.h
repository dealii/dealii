
#ifndef BROWSELINEEDIT_H
#define BROWSELINEEDIT_H

#include <QFrame>
#include <QLineEdit>
#include <QFileDialog>
#include <QPushButton>

/**
 * The BrowseLineEdit class provides a special line editor for the parameterGUI.
 * While editing file- or directory names it is much more easier to have a file-dialog
 * and just click on existing files or directories. This editor provides a simple QLineEditor
 * and a browse-button which opens a file- or a directory dialog. Clicking on existing files or directories
 * copies the path to the line editor. Depending on the <tt>BrowseType</tt> given in the constructor
 * the browse button opens a <tt>file</tt> or a <tt>directory</tt> dialog.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
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
				      * A <tt>slot</tt> for setting the text of the line editor.
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

#endif
