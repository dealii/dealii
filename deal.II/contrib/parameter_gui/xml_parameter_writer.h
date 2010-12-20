
#ifndef XMLPARAMETERWRITER_H
#define XMLPARAMETERWRITER_H

#include <QXmlStreamWriter>
#include <QTreeWidget>
#include <QTreeWidgetItem>


/**
 * The XMLParameterWriter class provides an interface to write parameters stored in a QTreeWidget to a file in XML format.
 * This class makes extensive use of the QXmlStreamWriter class, which implements the basic functionalities for writing
 * XML files.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
 */
class XMLParameterWriter
{
  public:
				     /**
				      * Constructor.
				      */
    XMLParameterWriter (QTreeWidget *tree_widget);
				     /**
				      * This function writes the parameters stored in the <tt>tree_widget</tt>
				      * to <tt>device</tt> in XML format. We use the QXmlStreaWriter class
				      * for this. The root element is
				      * <code><ParameterHandler></code>
				      */
    bool write_xml_file (QIODevice *device);

  private:
				     /**
				      * This function writes an <tt>item</tt> of <tt>tree_widget</tt>
				      * to a file in XML format. For this the QXmlStreamWriter class is used.
				      * If the <tt>item</tt> is a parameter, the elements that describes this parameter
				      * are written:
				      * @code
				      *   <value>value</value>
				      *   <default_value>default_value</default_value>
				      *   <documentation>documentation</documentation>
				      *   <pattern>pattern</pattern>
				      *   <pattern_description>[pattern_description]</pattern_description>
				      * @endcode
				      * If the <tt>item</tt> is a subsection, a start element <code>this_subsection</code> is written
				      * and <tt>write_item</tt> is called recursively to write the next <tt>item</tt>.
				      */
    void write_item (QTreeWidgetItem *item);
				     /**
				      * Reimplemented from the ParameterHandler class.
				      * Mangle a string so that it
				      * doesn't contain any special
				      * characters or spaces.
				      */
    QString  mangle (const QString &s);
				     /**
				      * An QXmlStreamWriter object
				      * which implements the functionalities
				      * we need for writing parameters to XML files.
				      */
    QXmlStreamWriter  xml;
				     /**
				      * A pointer to the QTreeWidget structure
				      * which stores the parameters.
				      */
    QTreeWidget * tree_widget;
};


#endif
