
#ifndef INFOMESSAGE_H
#define INFOMESSAGE_H

#include <QDialog>
#include <QSettings>
#include <QCheckBox>
#include <QTextEdit>
#include <QLabel>


/**
 * The InfoMessage class implements a special info message for the parameterGUI.
 * Besides showing a info message itself, the dialog shows a checkbox "Show this message again".
 * If the user unchecks this box, this is stored in the "settings.ini" file and will be reloaded
 * the next time the user opens the parameterGUI. The intention of such a info message is the following.
 * The user should have some information on how using the GUI "at hand"
 * such as "how to edit parameter values?" for example. But after reading this message, the user knows it
 * and the message should not appear permanently.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
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
				      * With this function the <tt>message</tt> which will be shown in the
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
    QAction *ok_action;
#endif
				     /**
				      * An object for storing <tt>settings</tt> in a file.
				      */
    QSettings * settings;
};

#endif
